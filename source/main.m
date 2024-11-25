function main(folderName, fileName, plotCheck, parFlag)
% folderName: location of input file, where output is placed as well
% fileName: name of input file, located in <folderName>
% plotCheck: must be set to 1, if plots are to be produced
% parFlag: 1 to use parfor loop for pencil decomposition

%% Simulation initialization and setup
% Parallelize with dual threaded cores (setup)
if parFlag
    % Separate parallel work folder prevents Matlab conflict on cluster,
    % when launching several cases in parallel.
    if ~exist([folderName 'matlabWork/'], 'dir')
      mkdir([folderName 'matlabWork/']);
    else
      rmdir([folderName 'matlabWork/'], 's');
      mkdir([folderName 'matlabWork/']);
    end
     c = parcluster;
%     c.NumWorkers = 3;
     c.JobStorageLocation = [folderName 'matlabWork/'];
     p = gcp('nocreate');
     if isempty(p)
        parpool(c);
     end
end

% Parse input files
fullName = [folderName, fileName];
[layerInfo, uniqueSpecies, uniqueSpeciesLatex, rxnInfo, faradaicRxnInfo, gasLayerInfo] = parseInputFile(fullName);
nSpecies = length(uniqueSpecies);
nGasSpecies = length(gasLayerInfo.gasSpeciesNames);

% Generate mesh - Uses generation function from multi-layer code.
[xCenter, xFace, dxC, dxF, yCenter, yFace, dyC, dyF, zCenter, zFace, dzC, dzF] = ...
    genOverallMesh(layerInfo(1).xInterfaces, ...
    [layerInfo(:).dxMax], [layerInfo(:).dxMin], ...
    [layerInfo(:).LID], [layerInfo(:).gridSymmetry], layerInfo(1).yInterfaces, layerInfo(1).dyMin, layerInfo(1).dyMax, layerInfo(1).yLRC, ...
    layerInfo(1).zInterfaces, layerInfo(1).dzMin, layerInfo(1).dzMax, layerInfo(1).zLRC);

% Gas modifications to remove ghost cell in x direction
gasLayerInfo.nFingers = floor(gasLayerInfo.spanwiseLength * gasLayerInfo.gasFingerRatio / (xFace(end-1) - xFace(2)));
gasLayerInfo.gasXYArea = gasLayerInfo.spanwiseLength ...
        * ((gasLayerInfo.gasDiffusionHeight + ((1-gasLayerInfo.gasFingerRatio) + gasLayerInfo.nFingers*(dxC(1)+dxC(end))/gasLayerInfo.spanwiseLength)*(yFace(end) - yFace(1)))*gasLayerInfo.gasPorosity ...
        + gasLayerInfo.gasChannelRatio*gasLayerInfo.gasChannelHeight);

% Physical constant array generation
constants = genPhysConstArrays(uniqueSpecies, layerInfo, rxnInfo, xCenter);

% Save some info for convenient post-processing
nX = length(xCenter);
nY = length(yCenter);
nZ = length(zCenter);
save([folderName, 'constants.mat'], 'constants');
save([folderName, 'uniqueSpecies.mat'], 'uniqueSpecies', 'uniqueSpeciesLatex');
save([folderName, 'meshSize.mat'], 'nX', 'nY', 'nZ', 'xFace', 'yFace', 'zFace', 'xCenter', 'yCenter', 'zCenter');
save([folderName, 'gasLayerInfo.mat'], 'gasLayerInfo');
save([folderName, 'rxnInfo.mat'], 'rxnInfo');
save([folderName, 'faradaicRxnInfo.mat'], 'faradaicRxnInfo');
save([folderName, 'layerInfo.mat'], 'layerInfo');

% Time paramters
% timeStepOverride = 1: uses parameters from input file
% timeStepOverride = 0: uses arbitrary parameters, listed below
if layerInfo(1).timeStepOverride
    tEnd = layerInfo(1).tEnd;
    dt = layerInfo(1).dt;
    dtOut = layerInfo(1).dtOut;
    dtPlot = layerInfo(1).dtPlot;
    tStart = layerInfo(1).tStart;
    
    if dt > min(layerInfo(1).dyMin, layerInfo(1).dzMin)^2/(max(constants.diff(1,:)) * 6)
        disp('Warning: dt is larger than explicit stability limit in y or z')
    end
else
    tEnd = 1e6;
    dt = 5e-1;
    dtOut = 1e4;
    dtPlot = 1e3;
    tStart = 0.00;
end
totalSteps = ceil((tEnd - tStart)/dt);
t = tStart;

% Variable formation and initialization
% Field variables are interlaced in the x-direction (first coordinate)
fields3D = zeros((nSpecies+1)*length(xCenter), length(yCenter), length(zCenter));

% Initialize gas field variable (1D, multiple species, eqipotential).
% Initialize using initial mole fraction, pressure throughout domain
gasFields = zeros(nZ * nGasSpecies, 1);
for sIndex = 1:nGasSpecies
    gasFields(sIndex:nGasSpecies:end) = gasLayerInfo.gasInitialMoleFraction(sIndex);
end
gasVelocity = ones(nZ+1, 1) * gasLayerInfo.gasFlowInlet / (gasLayerInfo.gasXYArea);

if dt > min(dzC)/max(gasVelocity)
    disp('Warning: dt is larger than explicit stability limit for velocity')
end

%% Optional homogeneous pre-equilibration of domain
% Homogeneous reaction equilibration (for outlet BC and initialization)
% Find equilibrium state, to initialize entire system and set the correct
% outflow/inflow Dirichlet concentration values.
if layerInfo(1).homogeneousEq
    stepCtr = 1;
    solveCounter = 0;
    tEnd_rxn = layerInfo(1).heq_tEnd;
    dt_rxn = layerInfo(1).heq_dt;
    tStart_rxn = 0;
    t_rxn = 0;
    totalStepsRxn = ceil((tEnd_rxn - tStart_rxn)/dt_rxn);

    % History term initialization
    dtEffOld = dt_rxn;
    dtEffOldOld = dtEffOld;

    eqConcValues = constants.initVal(1,:)';
    eqConcValuesOld = eqConcValues;
    eqConcValuesOldOld = eqConcValuesOld;
    eqConcValuesOldOldOld = eqConcValuesOldOld;

    for ii = (ceil(t_rxn/dt_rxn)+1):totalStepsRxn

        notConverged = 1;

        dtEff = dt_rxn;

        targetTime = tStart_rxn + ii*dt_rxn;

        while (notConverged) || abs(t_rxn - targetTime) > 10 * eps

            tEndEff = t_rxn + dtEff;
            eqConcValuesStar = eqConcValues;

            disp(' ')
            disp('###################################################')
            disp(['HOMOGENEOUS PRE-EQUILIBRATION. Step: ' num2str(solveCounter+1)])
            disp('###################################################')
            disp(' ')

            if ii == 1
                homogeneousRxnSourceTerms_init = computeRHS_makrand_homogeneousPre(nSpecies, eqConcValues, constants);
            end

            % Update concentration fields without homogeneous reactions
            concNorm = inf;
            oldConcNorm = -inf;
            concIterationCtr = 0;
            isImproving = 1;
            isFirstStep = 1;
            concNormInitial = -inf;

            while (concNorm > concNormInitial * 1e-9 && concIterationCtr < 10 && (abs(concNorm - oldConcNorm) > concNormInitial * 1e-12) && isImproving) || concIterationCtr < 3
                % Update concentration fields without homogeneous reactions
                homogeneousRxnSourceTerms = computeRHS_makrand_homogeneousPre(nSpecies, eqConcValuesStar, constants);

                [deltaVals, ~, residual] = doTimeStep_concentrationHomogeneousPre(eqConcValuesStar, eqConcValuesOld, eqConcValuesOldOld, eqConcValuesOldOldOld,...
                    t, dtEff, dtEffOld, dtEffOldOld, constants, homogeneousRxnSourceTerms, homogeneousRxnSourceTerms_init, stepCtr);

                eqConcValuesStar = eqConcValuesStar + deltaVals;

                oldConcNorm = concNorm;
                concNorm = sqrt(sum(residual.^2));

                disp(['Homogeneous pre-solve conc. |residual|: ' num2str(concNorm)]);

                if isFirstStep
                    concNormInitial = concNorm;
                    isFirstStep = 0;
                end

                if concNorm >= oldConcNorm
                    isImproving = 0;
                end
                concIterationCtr = concIterationCtr + 1;

            end

            % check whether any negative, nan, or imaginary values exist
            negVals = any(eqConcValuesStar < -eps);
            nanVals = any(isnan(eqConcValuesStar));
            imagVals = any(~isreal(eqConcValuesStar));

            % Check for charge accumulation
            netCharge = sum(constants.vale(1, :) .* (eqConcValuesStar - eqConcValuesOld)') / constants.m3ToLiter;
            bulkCharge = netCharge > layerInfo(1).chargeThreshold;

            % Force refinement flag, for debugging
            forceRefinement = 0;

            if nanVals || imagVals || forceRefinement %|| negVals  || bulkCharge
                if negVals
                    disp('REFINING DUE TO NEGATIVES')
                end
                if nanVals
                    disp('REFINING DUE TO NANS')
                end
                if imagVals
                    disp('REFINING DUE TO IMAGINARIES')
                end
                if bulkCharge
                    disp('REFINING DUE TO CHARGE')
                end
                notConverged = 1;
                dtEff = min([dtEff/2, 2*dtEffOld]);
                disp(['Smallest time step: ' num2str(dtEff)]);

            else
                t_rxn = tEndEff;

                notConverged = 0;

                solveCounter = solveCounter + 1;

                eqConcValues = eqConcValuesStar;
                eqConcValuesOldOldOld = eqConcValuesOldOld;
                eqConcValuesOldOld = eqConcValuesOld;
                eqConcValuesOld = eqConcValues;

                dtEffOldOld = dtEffOld;
                dtEffOld = dtEff;
                dtEff = min([targetTime - t_rxn, 1.1*dtEffOld]);

                disp(['Reached time: ' num2str(t_rxn)])
                stepCtr = stepCtr + 1;
            end
        end
    end

    constants.initVal = repmat(eqConcValues', [nX, 1]);
    save([folderName, 'constants.mat'], 'constants');
end

if layerInfo(1).yRightElectrolyteOverride
    layerInfo(1).yRightBCValues = constants.initVal(1,:);
end

%% Domain initialization
% Default behavior - initializes 1D profile, homogeneous in y and z
for yIndex = 1:size(fields3D, 2)
    for zIndex = 1:size(fields3D, 3)
        for ii = 1:nSpecies
            if layerInfo(1).yLeftNoFluxCondition == 0 && layerInfo(1).yRightNoFluxCondition == 0
                fields3D(ii:nSpecies+1:end, yIndex, zIndex) = layerInfo(1).yLeftBCValues(ii) + yCenter(yIndex)/yCenter(end)*(layerInfo(1).yRightBCValues(ii) - layerInfo(1).yLeftBCValues(ii));
            elseif layerInfo(1).zBottomNoFluxCondition == 0 && layerInfo(1).zTopNoFluxCondition == 0
                fields3D(ii:nSpecies+1:end, yIndex, zIndex) = layerInfo(1).zBottomBCValues(ii) + zCenter(zIndex)/zCenter(end)*(layerInfo(1).zTopBCValues(ii) - layerInfo(1).zBottomBCValues(ii));
            else
                if layerInfo(1).yLeftNoFluxCondition == 0
                    fields3D(ii:nSpecies+1:end, yIndex, zIndex) = layerInfo(1).yLeftBCValues(ii);
                elseif layerInfo(1).yRightNoFluxCondition == 0
                    fields3D(ii:nSpecies+1:end, yIndex, zIndex) = layerInfo(1).yRightBCValues(ii);
                elseif layerInfo(1).zBottomNoFluxCondition == 0
                    fields3D(ii:nSpecies+1:end, yIndex, zIndex) = layerInfo(1).zBottomBCValues(ii);
                elseif layerInfo(1).zTopNoFluxCondition == 0
                    fields3D(ii:nSpecies+1:end, yIndex, zIndex) = layerInfo(1).zTopBCValues(ii);
                else
                    fields3D(ii:(nSpecies+1):end, yIndex, zIndex) = constants.initVal(:,ii);
                    fields3D(ii, yIndex, zIndex) = layerInfo(1).xFrontBCValues(ii);
                    fields3D(end-(nSpecies+1)+ii, yIndex, zIndex) = layerInfo(1).xBackBCValues(ii);
                end
            end
        end
        
        if layerInfo(1).xFrontNoFluxCondition
            if layerInfo(end).xBackNoFluxCondition
                if layerInfo(1).yLeftNoFluxCondition == 0 && layerInfo(1).yRightNoFluxCondition == 0
                    fields3D(nSpecies+1:nSpecies+1:end, yIndex, zIndex) = layerInfo(1).yLeftPotentialBCValue + yCenter(yIndex)/yCenter(end)*(layerInfo(1).yRightPotentialBCValue - layerInfo(1).yLeftPotentialBCValue);
                elseif layerInfo(1).zBottomNoFluxCondition == 0 && layerInfo(1).zTopNoFluxCondition == 0
                    fields3D(nSpecies+1:nSpecies+1:end, yIndex, zIndex) = layerInfo(1).zBottomPotentialBCValue + zCenter(zIndex)/zCenter(end)*(layerInfo(1).zTopPotentialBCValue - layerInfo(1).zBottomPotentialBCValue);
                else
                    if layerInfo(1).yLeftNoFluxCondition == 0
                        fields3D(nSpecies+1:nSpecies+1:end, yIndex, zIndex) = layerInfo(1).yLeftPotentialBCValue;
                    elseif layerInfo(1).yRightNoFluxCondition == 0 || layerInfo(1).yRightElectrolyteOverride
                        fields3D(nSpecies+1:nSpecies+1:end, yIndex, zIndex) = layerInfo(1).yRightPotentialBCValue;
                    elseif layerInfo(1).zBottomNoFluxCondition == 0
                        fields3D(nSpecies+1:nSpecies+1:end, yIndex, zIndex) = layerInfo(1).zBottomPotentialBCValue;
                    elseif layerInfo(1).zTopNoFluxCondition == 0
                        fields3D(nSpecies+1:nSpecies+1:end, yIndex, zIndex) = layerInfo(1).zTopPotentialBCValue;
                    else
                        fields3D(nSpecies+1:nSpecies+1:end, yIndex, zIndex) = 0;
                    end
                end
            else
                fields3D(nSpecies+1:nSpecies+1:end, yIndex, zIndex) = ...
                    layerInfo(1).xBackPotentialBCValue;
            end
        else
            if layerInfo(end).xBackNoFluxCondition
                fields3D(nSpecies+1:nSpecies+1:end, yIndex, zIndex) = ...
                    layerInfo(1).xFrontPotentialBCValue;
            else
                fields3D(nSpecies+1:nSpecies+1:end, yIndex, zIndex) = ...
                    layerInfo(1).xFrontPotentialBCValue + xCenter/xCenter(end)...
                    * (layerInfo(1).xBackPotentialBCValue - layerInfo(1).xFrontPotentialBCValue);
                fields3D(nSpecies+1, yIndex, zIndex) = layerInfo(1).xFrontPotentialBCValue;
                fields3D(end, yIndex, zIndex) = layerInfo(1).xBackPotentialBCValue;
            end
        end
    end
end

% Create 3Dd-slice-form of BC vector for Dirichlet BC, in y direction
yLeftBCVec3DForm = zeros(size(fields3D(:, 1, :)));
for speciesIndex = 1:nSpecies
    for zIndex = 1:nZ
        yLeftBCVec3DForm(speciesIndex:(nSpecies+1):end, 1, zIndex) = layerInfo(1).yLeftBCValues(speciesIndex);
    end
end
% Neumann BC for potential
for zIndex = 1:nZ
    yLeftBCVec3DForm((nSpecies+1):(nSpecies+1):end, 1, zIndex) = layerInfo(1).yLeftPotentialBCValue;
end

yRightBCVec3DForm = zeros(size(fields3D(:, 1, :)));
for speciesIndex = 1:nSpecies
    for zIndex = 1:nZ
        yRightBCVec3DForm(speciesIndex:(nSpecies+1):end, 1, zIndex) = layerInfo(1).yRightBCValues(speciesIndex);
    end
end
% Neumann BC for potential
yRightBCVec3DForm((nSpecies+1):(nSpecies+1):end, 1, :) = layerInfo(1).yRightPotentialBCValue;
bcPotentialOld = layerInfo(1).yRightPotentialBCValue;

% Create 3Dd-slice-form of BC vector for Dirichlet BC, in z direction
zBottomBCVec3DForm = zeros(size(fields3D(:, :, 1)));
for speciesIndex = 1:nSpecies
    for yIndex = 1:nY
        zBottomBCVec3DForm(speciesIndex:(nSpecies+1):end, yIndex, 1) = layerInfo(1).zBottomBCValues(speciesIndex);
    end
end
% Neumann BC for potential
zBottomBCVec3DForm((nSpecies+1):(nSpecies+1):end, :, 1) = layerInfo(1).zBottomPotentialBCValue;

zTopBCVec3DForm = zeros(size(fields3D(:, :, 1)));
for speciesIndex = 1:nSpecies
    for yIndex = 1:nY
        zTopBCVec3DForm(speciesIndex:(nSpecies+1):end, yIndex, 1) = layerInfo(1).zTopBCValues(speciesIndex);
    end
end
% Neumann BC for potential
zTopBCVec3DForm((nSpecies+1):(nSpecies+1):end, :, 1) = layerInfo(1).zTopPotentialBCValue;

% Restart from old files
if ~layerInfo(1).newStart
    [fields3D, ~, ~, ~] = ...
        readBinaryOutput([folderName 'time_' num2str(layerInfo(1).restartFileTime) '.bin'], nX, nY, nZ, length(uniqueSpecies));
    
    fileID = fopen([folderName 'gasFields_time_' num2str(layerInfo(1).restartFileTime) '.bin']);
    gasFields = fread(fileID, 'double');
    fclose (fileID);
    
    fileID = fopen([folderName 'gasVelocity_time_' num2str(layerInfo(1).restartFileTime) '.bin']);
    gasVelocity = fread(fileID, 'double');
    fclose (fileID);
end

% Set gas inlet velocity (in case of restart with different parameter)
gasVelocity(1) = gasLayerInfo.gasFlowInlet / (gasLayerInfo.gasXYArea);

% Output frequency
outputInt = max(floor(dtOut/dt), 1);
plotInt = max(floor(dtPlot/dt), 1);

doPlot = plotCheck;
doOutputFunction = 1;       % 1 is on, will produce binary files
showUnconvergedPlots = 0;   % Show plots when time refinement is required

if doOutputFunction
    doOutput(xCenter, yCenter, zCenter, uniqueSpecies, fields3D, gasFields, gasVelocity, gasLayerInfo.gasSpeciesNames, t, folderName, bcPotentialOld);
end

% Plotting function here
if doPlot
    doAllPlots(fields3D((nSpecies+1)+1:end-(nSpecies+1), :, :), gasFields, gasVelocity, xCenter(2:end-1), yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, uniqueSpeciesLatex, gasLayerInfo.gasSpeciesNames, gasLayerInfo.gasSpeciesLatex, constants, gasLayerInfo)
end

% Initialize matrix for velocity solve
VelLHSOperator = [[1 zeros(1, nZ-1); diag(1./dzC(1:end-1), 1) - diag(1./dzC, 0)], [zeros(nZ,1); 1/dzC(end)]];

% Save original exchange current density, for linear ramp (to mitigate transient)
originalExchangeCurrents = zeros(length(faradaicRxnInfo), 1);
for rxnIndex = 1:length(faradaicRxnInfo)
    originalExchangeCurrents(rxnIndex) = faradaicRxnInfo(rxnIndex).exchangeCurrentDensity;
end

% Structure with cell volumes, for volume averaging
dvC= zeros(nX, nY, nZ);
for xIndex =1:nX
    for yIndex = 1:nY
        for zIndex = 1:nZ
            dvC(xIndex, yIndex, zIndex) = dxC(xIndex) * dyC(yIndex) * dzC(zIndex);
        end
    end
end

% Convenience broadcast variable initialization
totalCellsYZPlane = length(yCenter) * length(zCenter);
deltaFieldsForAccumulation = cell(totalCellsYZPlane, 1);
ohIndex = find(matches(uniqueSpecies, "OH-"));
hIndex = find(matches(uniqueSpecies, "H+"));
waterRxnIndices = layerInfo(1).waterReactionIndices;

if parFlag
    xStripMiddle_b = cell(totalCellsYZPlane, 1);
    xStripOld_b = cell(totalCellsYZPlane, 1);
    xStripOldOld_b = cell(totalCellsYZPlane, 1);
    xStripOldOldOld_b = cell(totalCellsYZPlane, 1);

    implicitTransportRHS_b = cell(totalCellsYZPlane, 1);
    homogeneousReactionRHS_b = cell(totalCellsYZPlane, 1);
    faradaicReactionRHS_b = cell(totalCellsYZPlane, 1);
    implicitTransportRHSOld_b = cell(totalCellsYZPlane, 1);
    explicitTransportRHSOld_b = cell(totalCellsYZPlane, 1);
    homogeneousReactionRHSOld_b = cell(totalCellsYZPlane, 1);
    faradaicReactionRHSOld_b = cell(totalCellsYZPlane, 1);
    explicitTransportRHSOldOld_b = cell(totalCellsYZPlane, 1);
    explicitTransportRHSOldOldOld_b = cell(totalCellsYZPlane, 1);

    xFrontNoFluxCondition = layerInfo(1).xFrontNoFluxCondition;
    xBackNoFluxCondition = layerInfo(1).xBackNoFluxCondition;
    isNeutralGas = layerInfo(1).isNeutralGas;
end

if layerInfo(1).galvanostaticMode
    measuredCurrentOld = computeDeviceCurrentDensity(fields3D, constants, dxC, dyC, dyF, dzC, nSpecies, gasLayerInfo, hIndex, ohIndex, waterRxnIndices, bcPotentialOld);
end

%% Time stepping loop
solveCounter = 0;

% History term initialization
dtForNextStep = dt;
dtEffOld = dt;
dtEffOldOld = dtEffOld;

fields3DOld = fields3D;
fields3DOldOld = fields3DOld;
fields3DOldOldOld = fields3DOldOld;

stepCtr = 1;

% Time stepping loop
for ii = (ceil(t/dt)+1):totalSteps
    
    notConverged = 1;
    if dtEffOldOld > dt / 3 || dtEffOld > dt / 2
        dtEff = dt; 
    else
        dtEff = dtForNextStep;
    end

    targetTime = tStart + ii*dt;
    adaptiveTimeCtr = 1;

    while (notConverged) || abs(t - targetTime) > 10 * eps
        tic
        tEndEff = t + dtEff;

        % Modify reaction rates according to linear ramp (avoids nasty
        % transient that produces super thin zones near boundary, due to rapid
        % equilibration of chemical reactions)
        for rxnIndex = 1:length(faradaicRxnInfo)
            faradaicRxnInfo(rxnIndex).exchangeCurrentDensity = min(tEndEff/layerInfo(1).reactionRampEndTime, 1) * originalExchangeCurrents(rxnIndex);
        end

        % If running in galvanostatic mode, adjust the potential boundary
        % condition accordingly, using target current and external
        % capacitance.
        if layerInfo(1).galvanostaticMode
            bcPotentialStar = bcPotentialOld - dtEff * (layerInfo(1).targetCurrent - measuredCurrentOld)/layerInfo(1).externalCapacitance;
            yRightBCVec3DForm(nSpecies+1:nSpecies+1:end, 1, :) = bcPotentialStar;

            if layerInfo(1).debugMode
                disp(['Measured current density: ' num2str(measuredCurrentOld)])
                disp(['Applied voltage: ' num2str(bcPotentialStar)])
            end
        end

        % "Initial" RHS terms are saved for Crank Nicolson startup
        if stepCtr == 1
            implicitTransportRHSOld = computeImplicitTransportRHS(nX, nY, nZ, dxC, dxF, nSpecies, fields3D, layerInfo(1).xFrontNoFluxCondition, layerInfo(1).xBackNoFluxCondition, constants, layerInfo(1).isNeutralGas);

            explicitTransportRHSOld = computeExplicitTransportRHS(dyC, dyF, dzC, dzF, nX, nY, nZ, fields3D, constants, ...
                layerInfo(1).yLeftNoFluxCondition, layerInfo(1).yRightNoFluxCondition, layerInfo(1).zBottomNoFluxCondition, layerInfo(1).zTopNoFluxCondition, ...
                yLeftBCVec3DForm, yRightBCVec3DForm, zBottomBCVec3DForm, zTopBCVec3DForm, layerInfo(1).yRightElectrolyteOverride, layerInfo(1).isNeutralGas, nSpecies, waterRxnIndices, ohIndex, hIndex);

            explicitTransportRHSOldOld = computeExplicitTransportRHS(dyC, dyF, dzC, dzF, nX, nY, nZ, fields3D, constants, ...
                layerInfo(1).yLeftNoFluxCondition, layerInfo(1).yRightNoFluxCondition, layerInfo(1).zBottomNoFluxCondition, layerInfo(1).zTopNoFluxCondition, ...
                yLeftBCVec3DForm, yRightBCVec3DForm, zBottomBCVec3DForm, zTopBCVec3DForm, layerInfo(1).yRightElectrolyteOverride, layerInfo(1).isNeutralGas, nSpecies, waterRxnIndices, ohIndex, hIndex);

            explicitTransportRHSOldOldOld = computeExplicitTransportRHS(dyC, dyF, dzC, dzF, nX, nY, nZ, fields3D, constants, ...
                layerInfo(1).yLeftNoFluxCondition, layerInfo(1).yRightNoFluxCondition, layerInfo(1).zBottomNoFluxCondition, layerInfo(1).zTopNoFluxCondition, ...
                yLeftBCVec3DForm, yRightBCVec3DForm, zBottomBCVec3DForm, zTopBCVec3DForm, layerInfo(1).yRightElectrolyteOverride, layerInfo(1).isNeutralGas, nSpecies, waterRxnIndices, ohIndex, hIndex);

            homogeneousReactionRHSOld = computeHomogeneousReactionRHS(nX, nY, nZ, nSpecies, fields3D, constants);

            faradaicReactionRHSOld = computeFaradaicReactionRHS(nX, nY, nZ, nSpecies, fields3D, faradaicRxnInfo, constants);
        end

        fields3DCoupled = fields3D;
        
        %%%%% SOLVE FOR GAS FIELD CONCENTRATIONS AND VELOCITY
        % Gas velocity ODE solve
        gasRHSCellCenters = zeros(nZ, 1);
        for sIndex = 1:nGasSpecies
            [speciesFaradaicReactionRate, speciesDissolutionRate] = computeGasProductionRate(fields3DCoupled, dxC, dxF, dyC, dyF, dzC, dzF, nX, nY, nZ, gasLayerInfo, constants, uniqueSpecies, faradaicRxnInfo, sIndex);
            gasRHSCellCenters = gasRHSCellCenters + (speciesFaradaicReactionRate + speciesDissolutionRate) * constants.NA * constants.kb * constants.T / (gasLayerInfo.gasPressure);
        end
        gasRHS = [gasVelocity(1); gasRHSCellCenters];
        
        gasVelocity = VelLHSOperator \ gasRHS;
        
        % RK4 method, explicit
        k1Frac = dtEff * rhsFunctionGasConc(gasFields, gasVelocity, fields3DCoupled, dxC, dxF, dyC, dyF, dzC, dzF, nX, nY, nZ, gasLayerInfo, constants, uniqueSpecies, faradaicRxnInfo);
        
        k2Frac = dtEff * rhsFunctionGasConc(gasFields + 1/2*k1Frac, gasVelocity, fields3DCoupled, dxC, dxF, dyC, dyF, dzC, dzF, nX, nY, nZ, gasLayerInfo, constants, uniqueSpecies, faradaicRxnInfo);
        
        k3Frac = dtEff * rhsFunctionGasConc(gasFields + 1/2*k2Frac, gasVelocity, fields3DCoupled, dxC, dxF, dyC, dyF, dzC, dzF, nX, nY, nZ, gasLayerInfo, constants, uniqueSpecies, faradaicRxnInfo);
        
        k4Frac = dtEff * rhsFunctionGasConc(gasFields + k3Frac, gasVelocity, fields3DCoupled, dxC, dxF, dyC, dyF, dzC, dzF, nX, nY, nZ, gasLayerInfo, constants, uniqueSpecies, faradaicRxnInfo);
        
        gasFieldsStar = gasFields ...
            + 1/6 * k1Frac ...
            + 1/3 * (k2Frac + k3Frac) ...
            + 1/6 * (k4Frac);
        
        % Henry's Law implementation for neutralGas species
        for sIndex = 1:nGasSpecies
            if layerInfo(1).isNeutralGas(sIndex)
                aqueousIndex = gasLayerInfo.isTrackedInElectrolyte(sIndex);
                for zIndex = 1:nZ
                    fields3DCoupled(aqueousIndex, :, zIndex) = (gasLayerInfo.henrysLawConstant(sIndex) * gasFieldsStar((zIndex-1)*nGasSpecies + sIndex) * gasLayerInfo.gasPressure);
                    fields3DCoupled(end-(nSpecies+1)+aqueousIndex, :, zIndex) = (gasLayerInfo.henrysLawConstant(sIndex) * gasFieldsStar((zIndex-1)*nGasSpecies + sIndex) * gasLayerInfo.gasPressure);
                    yLeftBCVec3DForm(aqueousIndex:(nSpecies+1):end, 1, zIndex) = (gasLayerInfo.henrysLawConstant(sIndex) * gasFieldsStar((zIndex-1)*nGasSpecies + sIndex) * gasLayerInfo.gasPressure);
                end
            end
        end

        deltaFields = zeros(size(fields3D));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%   Homogeneous and Faradaic Rxn.  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(' ')
        disp('######################################################')
        disp(['[HOMOGENEOUS EQUILIBRATION ITERATIONS. Step: ' num2str(solveCounter+1)])
        disp('#######################################################')
        disp(' ')

        blockNormInitial = 0;
        oldBlockNorm = -inf;
        blockNorm = inf;
        blockIterationCtr = 0;
        isFirstBlockIteration = 1;

        while (blockNorm > blockNormInitial * 1e-6 && blockIterationCtr < 10 && abs(blockNorm - oldBlockNorm) > blockNormInitial * 1e-2) || blockIterationCtr < 2
            % Preparing variables for concentration solve
            concNorm = inf;
            oldConcNorm = -inf;
            concNormInitial = 0;
            isFirstStep = 1;
            concIterationCtr = 0;
            isImproving = 1;
            
            % Concentration solve should converge in one step due to
            % linearity, need to debug and see why this is not happening.
            while (concNorm > concNormInitial * 1e-6 && concIterationCtr < 10 && (abs(concNorm - oldConcNorm) > concNormInitial * 1e-2) && isImproving) || concIterationCtr < 2
                % Update concentration fields without homogeneous reactions
                implicitTransportRHS = computeImplicitTransportRHS(nX, nY, nZ, dxC, dxF, nSpecies, fields3DCoupled, layerInfo(1).xFrontNoFluxCondition, layerInfo(1).xBackNoFluxCondition, constants, layerInfo(1).isNeutralGas);

                homogeneousReactionRHS = computeHomogeneousReactionRHS(nX, nY, nZ, nSpecies, fields3DCoupled, constants);

                faradaicReactionRHS = computeFaradaicReactionRHS(nX, nY, nZ, nSpecies, fields3DCoupled, faradaicRxnInfo, constants);
                
                %%%%%% PARALLELIZE HERE FOR SPEED %%%%%%
                concResidual = zeros(nY, nZ);
                if parFlag
                    for indexYZPlane = 1:totalCellsYZPlane
                        [yIndex, zIndex] = ind2sub([length(yCenter), length(zCenter)], indexYZPlane);
                        xStripMiddle_b{indexYZPlane} = fields3DCoupled(:, yIndex, zIndex);
                        xStripOld_b{indexYZPlane} = fields3DOld(:, yIndex, zIndex);
                        xStripOldOld_b{indexYZPlane} = fields3DOldOld(:, yIndex, zIndex);
                        xStripOldOldOld_b{indexYZPlane} = fields3DOldOldOld(:, yIndex, zIndex);
                        implicitTransportRHS_b{indexYZPlane} = implicitTransportRHS(:, yIndex, zIndex, :);
                        homogeneousReactionRHS_b{indexYZPlane} = homogeneousReactionRHS(:, yIndex, zIndex, :);
                        faradaicReactionRHS_b{indexYZPlane} = faradaicReactionRHS(:, yIndex, zIndex, :);
                        implicitTransportRHSOld_b{indexYZPlane} = implicitTransportRHSOld(:, yIndex, zIndex, :);
                        explicitTransportRHSOld_b{indexYZPlane} = explicitTransportRHSOld(:, yIndex, zIndex, :);
                        homogeneousReactionRHSOld_b{indexYZPlane} = homogeneousReactionRHSOld(:, yIndex, zIndex, :);
                        faradaicReactionRHSOld_b{indexYZPlane} = faradaicReactionRHSOld(:, yIndex, zIndex, :);
                        explicitTransportRHSOldOld_b{indexYZPlane} = explicitTransportRHSOldOld(:, yIndex, zIndex, :);
                        explicitTransportRHSOldOldOld_b{indexYZPlane} = explicitTransportRHSOldOldOld(:, yIndex, zIndex, :);
                    end

                    parfor indexYZPlane = 1:totalCellsYZPlane
                        [deltaXStripOut, ~, residual] = doTimeStep_concentrationWithHomogeneous(dxC, dxF, nSpecies, xStripMiddle_b{indexYZPlane}, xStripOld_b{indexYZPlane}, xStripOldOld_b{indexYZPlane}, xStripOldOldOld_b{indexYZPlane}, xFrontNoFluxCondition, xBackNoFluxCondition, ...
                            t, dtEff, dtEffOld, dtEffOldOld, constants, faradaicRxnInfo, isNeutralGas, ...
                            implicitTransportRHS_b{indexYZPlane}, homogeneousReactionRHS_b{indexYZPlane}, faradaicReactionRHS_b{indexYZPlane}, ...
                            implicitTransportRHSOld_b{indexYZPlane}, explicitTransportRHSOld_b{indexYZPlane}, homogeneousReactionRHSOld_b{indexYZPlane}, faradaicReactionRHSOld_b{indexYZPlane}, ...
                            explicitTransportRHSOldOld_b{indexYZPlane}, explicitTransportRHSOldOldOld_b{indexYZPlane}, stepCtr);
                        
                        deltaFieldsForAccumulation{indexYZPlane} = deltaXStripOut;
        
                        concResidual(indexYZPlane) = sum(residual.^2);
                    end
                else
                    for indexYZPlane = 1:totalCellsYZPlane
                        [yIndex, zIndex] = ind2sub([length(yCenter), length(zCenter)], indexYZPlane);
                        
                        xStripMiddle = fields3DCoupled(:, yIndex, zIndex);
                        xStripOld = fields3DOld(:, yIndex, zIndex);
                        xStripOldOld = fields3DOldOld(:, yIndex, zIndex);
                        xStripOldOldOld = fields3DOldOldOld(:, yIndex, zIndex);
                            
                        [deltaXStripOut, ~, residual] = doTimeStep_concentrationWithHomogeneous(dxC, dxF, nSpecies, xStripMiddle, xStripOld, xStripOldOld, xStripOldOldOld, layerInfo(1).xFrontNoFluxCondition, layerInfo(end).xBackNoFluxCondition, ...
                            t, dtEff, dtEffOld, dtEffOldOld, constants, faradaicRxnInfo, layerInfo(1).isNeutralGas, ...
                            implicitTransportRHS(:, yIndex, zIndex, :), homogeneousReactionRHS(:, yIndex, zIndex, :), faradaicReactionRHS(:, yIndex, zIndex, :), ...
                            implicitTransportRHSOld(:, yIndex, zIndex, :), explicitTransportRHSOld(:, yIndex, zIndex, :), homogeneousReactionRHSOld(:, yIndex, zIndex, :), faradaicReactionRHSOld(:, yIndex, zIndex, :), ...
                            explicitTransportRHSOldOld(:, yIndex, zIndex, :), explicitTransportRHSOldOldOld(:, yIndex, zIndex, :), stepCtr);
                        
                        deltaFieldsForAccumulation{indexYZPlane} = deltaXStripOut;
        
                        concResidual(indexYZPlane) = sum(residual.^2);
                    end
                end
                
                oldConcNorm = concNorm;
                concNorm = sqrt(sum(concResidual, [1, 2]));
                
                if layerInfo(1).debugMode
                    disp(['Homogeneous Iteration ' num2str(concIterationCtr) ' conc. |residual|: ' num2str(concNorm, '%.4e')]);
                end

                if isFirstStep
                    concNormInitial = concNorm;
                    oldBlockNorm = blockNorm;
                    blockNorm = concNorm;
                    isFirstStep = 0;
                end

                if isFirstBlockIteration
                    blockNormInitial = concNorm;
                    isFirstBlockIteration = 0;
                end

                % Update concentration fields
                for indexYZPlane = 1:totalCellsYZPlane
                    [yIndex, zIndex] = ind2sub([length(yCenter), length(zCenter)], indexYZPlane);
                    for sIndex = 1:nSpecies
                        deltaFields(sIndex:nSpecies+1:end, yIndex, zIndex) = ...
                            [0; deltaFieldsForAccumulation{indexYZPlane}(sIndex:nSpecies:end); 0];
                    end
                end
                fields3DCoupled = fields3DCoupled + deltaFields;

                if concNorm >= oldConcNorm
                    isImproving = 0;
                end
                
                concIterationCtr = concIterationCtr + 1;

            end           

            % Update electric potential (full nonlinear iteration)
            potentialNorm = inf;
            oldPotentialNorm = -inf;
            potentialNormInitial = 0;
            potentialIterationCtr = 0;
            potentialMaxDelta = inf;
            while potentialNorm > potentialNormInitial * 1e-6 && potentialIterationCtr < 10 && (abs(potentialNorm - oldPotentialNorm) > potentialNormInitial * 1e-2) && potentialMaxDelta > eps(1)
                % Update electric potential
                oldPotentialNorm = potentialNorm;
                
                [deltaPotential, potentialNorm] = doTimeStep_potentialCoupling(dxC, dxF, dyC, dyF, dzC, dzF, ...
                    nSpecies, fields3DCoupled,layerInfo(1).xFrontNoFluxCondition, layerInfo(end).xBackNoFluxCondition, layerInfo(1).yLeftNoFluxCondition, layerInfo(1).yRightNoFluxCondition, layerInfo(1).zBottomNoFluxCondition, layerInfo(1).zTopNoFluxCondition, ...
                    constants, yLeftBCVec3DForm, yRightBCVec3DForm, zBottomBCVec3DForm, zTopBCVec3DForm, faradaicRxnInfo, layerInfo(1).isNeutralGas, layerInfo(1).yRightElectrolyteOverride, waterRxnIndices, ohIndex, hIndex);

                fields3DCoupled(nSpecies+1:nSpecies+1:end, :, :) = ...
                    fields3DCoupled(nSpecies+1:nSpecies+1:end, :, :) + deltaPotential;

                if potentialIterationCtr == 0
                    potentialNormInitial = potentialNorm;
                end

                potentialMaxDelta = max(abs(deltaPotential(:)));
                if layerInfo(1).debugMode
                    disp(['Homogeneous Iteration ' num2str(potentialIterationCtr) ' potential |residual|: ' num2str(potentialNorm, '%.4e')]);
                end

                potentialIterationCtr = potentialIterationCtr + 1;
            end

            if layerInfo(1).debugMode   
                netCharge = zeros(nX-2, nY, nZ);
                for sIndex = 1:nSpecies
                    netCharge = netCharge + constants.vale(1, sIndex) * (fields3DCoupled(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), :, :) - fields3D(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), :, :)) / constants.m3ToLiter;
                end
                
                integratedCharge = sqrt(sum(netCharge .^2 .* dvC(2:end-1, :, :), [1 2 3]) / (xFace(end)-xFace(1)) / (yFace(end)-yFace(1)) / (zFace(end)-zFace(1)));
                
                disp(' ')
                disp(['Homogeneos + Faradaic iteration ' num2str(blockIterationCtr) ', Net charge generated: ' num2str(integratedCharge)])
                disp(' ')
            end

            blockIterationCtr = blockIterationCtr + 1;

            blockNorms = zeros(nSpecies+1, 1);
            for sIndex = 1:nSpecies+1
                blockNorms(sIndex) = norm(reshape(deltaFields(sIndex:nSpecies+1:end, :, :), [nX*nY*nZ, 1])) / norm(reshape(fields3DCoupled(sIndex:nSpecies+1:end, :, :), [nX*nY*nZ, 1]));
            end
            
            if layerInfo(1).debugMode
                disp(' ')
                disp(['Norm of block iteration delta: ' num2str(max(blockNorms), '%.4e')])
                disp(' ')
            end
        end

        fields3DStar = fields3DCoupled;
        
        % check whether any negative values exist
        negVals = 0;
        for kk = 1:nSpecies
            negVals = negVals || any(any(any(fields3DStar((nSpecies+1)+kk:(nSpecies+1):end-(nSpecies+1), :, :) <- eps)));
        end
        
        nanVals = any(any(any(isnan(fields3DStar))));
        imagVals = any(any(any(~isreal(fields3DStar))));

        netCharge = zeros(nX-2, nY, nZ);
        for sIndex = 1:nSpecies
            netCharge = netCharge + constants.vale(1, sIndex) * (fields3DStar(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), :, :) - fields3D(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), :, :)) / constants.m3ToLiter;
        end
        integratedCharge = sqrt(sum(netCharge .^2 .* dvC(2:end-1, :, :), [1 2 3]) / (xFace(end)-xFace(1)) / (yFace(end)-yFace(1)) / (zFace(end)-zFace(1)));
        
        bulkCharge = integratedCharge > layerInfo(1).chargeThreshold;

        forceRefinement = 0;
        
        if negVals || nanVals || imagVals || forceRefinement || bulkCharge
            if layerInfo(1).debugMode
                disp(' ')
                disp('****************************')
                if negVals
                    disp('REFINING DUE TO NEGATIVES')
                end
                if nanVals
                    disp('REFINING DUE TO NANS')
                end
                if imagVals
                    disp('REFINING DUE TO IMAGINARIES')
                end
                if bulkCharge
                    disp('REFINING DUE TO CHARGE')
                end
                disp('****************************')
                disp(' ')
            end
            
            notConverged = 1;
            
            if adaptiveTimeCtr == 1
                dtEff = dtEff/2;
            else
                dtEff = 0.9 * dtEff;
            end
            
            disp(['Refining time step to: ' num2str(dtEff)]);
            
            if doPlot && showUnconvergedPlots
                doAllPlots(fields3DStar, gasFieldsStar, gasVelocity, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, uniqueSpeciesLatex, gasLayerInfo.gasSpeciesNames, gasLayerInfo.gasSpeciesLatex, constants, gasLayerInfo)
            end
            
        else
            t = tEndEff;

            notConverged = 0;
            solveCounter = solveCounter + 1;

            gasFields = gasFieldsStar;
            
            fields3D = fields3DStar;
            fields3DOldOldOld = fields3DOldOld;
            fields3DOldOld = fields3DOld;
            fields3DOld = fields3D;

            dtEffOldOld = dtEffOld;
            dtEffOld = dtEff;
            if integratedCharge < 0.80*layerInfo(1).chargeThreshold
                dtEff = min([targetTime - t, 1.1*dtEffOld]);
                dtForNextStep = 1.1*dtEffOldOld;
            else
                dtEff = min([targetTime - t, dtEffOld]);
                dtForNextStep = dtEffOldOld;
            end

            netCharge = zeros(nX-2, nY, nZ);
            for sIndex = 1:nSpecies
                netCharge = netCharge + constants.vale(1, sIndex) * (fields3D(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), :, :)) / constants.m3ToLiter;
            end
            
            integratedCharge = sqrt(sum(netCharge .^ 2 .* dvC(2:end-1, :, :), [1 2 3]) / (xFace(end)-xFace(1)) / (yFace(end)-yFace(1)) / (zFace(end)-zFace(1)));

            implicitTransportRHSOld = computeImplicitTransportRHS(nX, nY, nZ, dxC, dxF, nSpecies, fields3D, layerInfo(1).xFrontNoFluxCondition, layerInfo(1).xBackNoFluxCondition, constants, layerInfo(1).isNeutralGas);
            explicitTransportRHSOldOldOld = explicitTransportRHSOldOld;
            explicitTransportRHSOldOld = explicitTransportRHSOld;
            explicitTransportRHSOld = computeExplicitTransportRHS(dyC, dyF, dzC, dzF, nX, nY, nZ, fields3D, constants, ...
                layerInfo(1).yLeftNoFluxCondition, layerInfo(1).yRightNoFluxCondition, layerInfo(1).zBottomNoFluxCondition, layerInfo(1).zTopNoFluxCondition, ...
                yLeftBCVec3DForm, yRightBCVec3DForm, zBottomBCVec3DForm, zTopBCVec3DForm, layerInfo(1).yRightElectrolyteOverride, layerInfo(1).isNeutralGas, nSpecies, waterRxnIndices, ohIndex, hIndex);
            homogeneousReactionRHSOld = computeHomogeneousReactionRHS(nX, nY, nZ, nSpecies, fields3D, constants);
            faradaicReactionRHSOld = computeFaradaicReactionRHS(nX, nY, nZ, nSpecies, fields3D, faradaicRxnInfo, constants);

            if layerInfo(1).galvanostaticMode
                bcPotentialOld = bcPotentialStar;
                measuredCurrentOld = computeDeviceCurrentDensity(fields3DOld, constants, dxC, dyC, dyF, dzC, nSpecies, gasLayerInfo, hIndex, ohIndex, waterRxnIndices, bcPotentialOld);
            end

            disp(' ')
            disp('---------------------------------------------------------')
            disp('---------------------------------------------------------')
            disp(['Reached time: ' num2str(t) ' using time step of ' num2str(dtEffOld) '.'])
            if layerInfo(1).debugMode
                disp(['RMS charge at end of time step: ' num2str(integratedCharge)])
            end
            disp('---------------------------------------------------------')
            disp('---------------------------------------------------------')
            disp(' ')

            stepCtr = stepCtr + 1;
            adaptiveTimeCtr = adaptiveTimeCtr + 1;

%             doAllPlots(fields3D, gasFields, gasVelocity, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, uniqueSpeciesLatex, gasLayerInfo.gasSpeciesNames, gasLayerInfo.gasSpeciesLatex, constants, gasLayerInfo)
        end
        toc
    end
    
    t = targetTime;
    
     if (mod(ii, outputInt) == 0) && (doOutputFunction)
        % output generation function here
        doOutput(xCenter, yCenter, zCenter, uniqueSpecies, fields3D, gasFields, gasVelocity, gasLayerInfo.gasSpeciesNames, t, folderName, bcPotentialOld);
        disp(['Completed step ' num2str(ii) ' of ' num2str(totalSteps) ' (' num2str(100*ii/totalSteps) '%)'])
    end
    
    if (mod(ii, plotInt) == 0) && doPlot
        % Plotting function here
        doAllPlots(fields3D, gasFields, gasVelocity, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, uniqueSpeciesLatex, gasLayerInfo.gasSpeciesNames, gasLayerInfo.gasSpeciesLatex, constants, gasLayerInfo)
    end
end
end
