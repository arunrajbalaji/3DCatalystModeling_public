function fullMMS(folderName, fileName, plotCheck)
% folderName: location of output
% fileName: name of input file, located in <folderName>
% plotCheck: must be set to 1, if plots are to be produced

% Parallelize with dual threaded cores
p = gcp('nocreate');

if ~exist([folderName 'matlabWork/'], 'dir') && isempty(p)
  mkdir([folderName 'matlabWork/']);
elseif isempty(p)
  rmdir([folderName 'matlabWork/'], 's');
  mkdir([folderName 'matlabWork/']);
end

% if isempty(p)
%     c = parcluster;
%     c.NumWorkers = 24;
%     c.JobStorageLocation = [folderName 'matlabWork/'];
%     parpool(c, 24);
% end

% addpath('~/multiLayerEPNP/advanpix/');
% mp.Digits(32);
% format longG;

% Parse input files
fullName = [folderName, fileName];

[layerInfo, uniqueSpecies, rxnInfo, faradaicRxnInfo, gasLayerInfo] = parseInputFile(fullName);

nSpecies = length(uniqueSpecies);
nGasSpecies = length(gasLayerInfo.gasSpeciesNames);

% Generate mesh - Need some way to determine dx limits, LRC
% Function for determining dx limits here
% LRC should be provided as part of input
[xCenter, xFace, dxC, dxF, yCenter, yFace, dyC, dyF, zCenter, zFace, dzC, dzF] = ...
    genOverallMesh(layerInfo(1).xInterfaces,...
    [layerInfo(:).dxMax], [layerInfo(:).dxMin],...cd highC
    [layerInfo(:).LID], [layerInfo(:).gridSymmetry], layerInfo(1).yInterfaces, layerInfo(1).dyMin, layerInfo(1).dyMax, layerInfo(1).yLRC, ...
    layerInfo(1).zInterfaces, layerInfo(1).dzMin, layerInfo(1).dzMax, layerInfo(1).zLRC);

% Gas modifications to remove ghost cell in x direction
gasLayerInfo.nFingers = floor(gasLayerInfo.spanwiseLength * gasLayerInfo.gasFingerRatio / (xFace(end-1) - xFace(2)));
gasLayerInfo.gasXYArea = gasLayerInfo.spanwiseLength ...
        * ((gasLayerInfo.gasDiffusionHeight + ((1-gasLayerInfo.gasFingerRatio) - gasLayerInfo.nFingers*(dxC(1)+dxC(end))/gasLayerInfo.spanwiseLength)*(yFace(end) - yFace(1)))*gasLayerInfo.gasPorosity ...
        + gasLayerInfo.gasChannelRatio*gasLayerInfo.gasChannelHeight);

% Mesh variables for plotting
[xCenterPlot, yCenterPlot, zCenterPlot] = meshgrid(xCenter, yCenter, zCenter);

% Physical constant array generation
constants = genPhysConstArrays(uniqueSpecies, layerInfo, rxnInfo, xCenter);

% Save some info for convenient post-processing
nX = length(xCenter);
nY = length(yCenter);
nZ = length(zCenter);
save([folderName, 'constants.mat'], 'constants');
save([folderName, 'uniqueSpecies.mat'], 'uniqueSpecies');
save([folderName, 'meshSize.mat'], 'nX', 'nY', 'nZ', 'xFace', 'yFace', 'zFace', 'xCenter', 'yCenter', 'zCenter');
save([folderName, 'gasLayerInfo.mat'], 'gasLayerInfo');
save([folderName, 'rxnInfo.mat'], 'rxnInfo');
save([folderName, 'faradaicRxnInfo.mat'], 'faradaicRxnInfo');
save([folderName, 'layerInfo.mat'], 'layerInfo');

% Time paramters, maybe start of simInfo eventually,
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

% Variable formation and initialization function here
% Field variables are interlaced in the x-direction (first coordinate)
fields3D = zeros((nSpecies+1)*length(xCenter), length(yCenter), length(zCenter));

% Initialize gas field variable (1D, multiple species, eqipotential).
% Initialize using initial mole fraction, pressure throughout domain
gasFields = zeros(nZ * nGasSpecies, 1);
for sIndex = 1:nGasSpecies
    gasFields(sIndex:nGasSpecies:end) = gasLayerInfo.gasInitialMoleFraction(sIndex);
    % gasFields(sIndex:nGasSpecies:end) = gasLayerInfo.gasInitialMoleFraction(sIndex) ...
    %         * gasLayerInfo.gasPressure/(constants.kb * constants.NA * constants.T) * cos(linspace(0, pi/2, nZ));
end
gasVelocity = ones(nZ+1, 1) * gasLayerInfo.gasFlowInlet / ...
    (gasLayerInfo.gasXYArea) / gasLayerInfo.gasTortuosity;

if dt > min(dzC)/max(gasVelocity)
    disp('Warning: dt is larger than explicit stability limit for velocity')
end

% Outer bulk electrolyte parameter calcualation
% bulkElectrolyteL = layerInfo(1).electrolyteKinematicViscosity^(1/6) * gasLayerInfo.spanwiseLength^(1/2) ...
%     * ((constants.diff(1,1) + constants.diff(1,4))/2)^(1/3) ...
%     / (layerInfo(1).electrolyteFlowRate / gasLayerInfo.spanwiseLength / layerInfo(1).bulkElectrolyteThickness)^(1/2);
bulkElectrolyteL = 1e10;

% Reading voltage sweep function if necessary
if layerInfo(1).voltageSweepOn
    voltageVals = layerInfo(1).voltageValues;
    voltageRate = layerInfo(1).voltageSweepRate;
    timeVals = tStart:(1/voltageRate):(length(voltageVals)-1)/voltageRate;
    
    layerInfo(1).yRightPotentialBCValue = voltageVals(1);
end

% Homogeneous reaction equilibration (for outlet BC and initialization)
% Find equilibrium state, to initialize entire system and set the correct
% outflow/inflow Dirichlet concentration values.

% Iteration variables
stepCtr = 1;
solveCounter = 0;

% Time stepping parameters and duration parameters
solveCounter = 0;
tEnd_rxn = 20;
dt_rxn = 1;
tStart_rxn = 0;
t_rxn = 0;
totalStepsRxn = ceil((tEnd_rxn - tStart_rxn)/dt_rxn);

% History term initialization
dtEffOld = dt_rxn;
dtEffOldOld = dtEffOld;

eqConcValues = constants.initVal(1,:)' * constants.m3ToLiter;
eqConcValuesOld = eqConcValues;
eqConcValuesOldOld = eqConcValuesOld;
eqConcValuesOldOldOld = eqConcValuesOldOld;

rxns = constants.rxns{1};
mmsEqVals = [7e-3; 1; 0; 0; 0; 3; 1] * constants.m3ToLiter;
mmsEpsilon = 1e-2;
mmsPeriod = 1;
mmsConc = @(mmsT, sIndex) mmsEqVals(sIndex) * (1 + mmsEpsilon * sin(2*pi*mmsT / mmsPeriod));
mmsTimeTerm = @(mmsT) constants.poro(1) * mmsEqVals * mmsEpsilon * 2 * pi / mmsPeriod * cos(2 * pi * mmsT / mmsPeriod);
mmsReactionTerm = @(mmsT) ...
    -constants.poro(1) * [-rxns(5, 5) * mmsConc(mmsT, 1) * mmsConc(mmsT, 7) + rxns(6, 5) * mmsConc(mmsT, 5) - rxns(7, 5) * mmsConc(mmsT, 1) + rxns(8, 5) * mmsConc(mmsT, 4); ...
    rxns(3, 5) * mmsConc(mmsT, 5) - rxns(4, 5) * mmsConc(mmsT, 2) * mmsConc(mmsT, 3) + rxns(9, 5) * mmsConc(mmsT, 5) * mmsConc(mmsT, 5) - rxns(10, 5) * mmsConc(mmsT, 2) * mmsConc(mmsT, 4) + rxns(13, 5) * mmsConc(mmsT, 5) * mmsConc(mmsT, 7) - rxns(14, 5) * mmsConc(mmsT, 2); ...
    rxns(1, 5) * mmsConc(mmsT, 4) - rxns(2, 5) * mmsConc(mmsT, 3) * mmsConc(mmsT, 5) + rxns(3, 5) * mmsConc(mmsT, 5) - rxns(4, 5) * mmsConc(mmsT, 2) * mmsConc(mmsT, 3) + rxns(15, 5) - rxns(16, 5) * mmsConc(mmsT, 3) * mmsConc(mmsT, 7); ...
    -rxns(1, 5) * mmsConc(mmsT, 4) + rxns(2, 5) * mmsConc(mmsT, 5) * mmsConc(mmsT, 3) + rxns(7, 5) * mmsConc(mmsT, 1) - rxns(8, 5) * mmsConc(mmsT, 4) + rxns(9, 5) * mmsConc(mmsT, 5) * mmsConc(mmsT, 5) - rxns(10, 5) * mmsConc(mmsT, 2) * mmsConc(mmsT, 4) - rxns(11, 5) * mmsConc(mmsT, 4) * mmsConc(mmsT, 7) + rxns(12, 5) * mmsConc(mmsT, 5); ...
    rxns(1, 5) * mmsConc(mmsT, 4) - rxns(2, 5) * mmsConc(mmsT, 5) * mmsConc(mmsT, 3) - rxns(3, 5) * mmsConc(mmsT, 5) + rxns(4, 5) * mmsConc(mmsT, 2) * mmsConc(mmsT, 3) + rxns(5, 5) * mmsConc(mmsT, 1) * mmsConc(mmsT, 7) - rxns(6, 5) * mmsConc(mmsT, 5) - rxns(9, 5) * mmsConc(mmsT, 5) * mmsConc(mmsT, 5) + rxns(10, 5) * mmsConc(mmsT, 2) * mmsConc(mmsT, 4) + rxns(11, 5) * mmsConc(mmsT, 4) * mmsConc(mmsT, 7) - rxns(12, 5) * mmsConc(mmsT, 5) - rxns(13, 5) * mmsConc(mmsT, 5) * mmsConc(mmsT, 7) + rxns(14, 5) * mmsConc(mmsT, 2); ...
    0; ...
    -rxns(5, 5) * mmsConc(mmsT, 1) * mmsConc(mmsT, 7) + rxns(6, 5) * mmsConc(mmsT, 5) - rxns(11, 5) * mmsConc(mmsT, 4) * mmsConc(mmsT, 7) + rxns(12, 5) * mmsConc(mmsT, 5) - rxns(13, 5) * mmsConc(mmsT, 5) * mmsConc(mmsT, 7) + rxns(14, 5) * mmsConc(mmsT, 2) + rxns(15, 5) - rxns(16, 5) * mmsConc(mmsT, 3) * mmsConc(mmsT, 7)];

for ii = (ceil(t_rxn/dt_rxn)+1):totalStepsRxn
    
    notConverged = 1;
    
    dtEff = dt_rxn;
%     if ii == 1
%         dtEff = 1e-11;
%     end
    
    targetTime = t_rxn+dt_rxn;
    
    while (notConverged) || (t_rxn ~= targetTime)

        tEndEff = t_rxn + dtEff;
        eqConcValuesStar = eqConcValues;

        mmsForcingTermTotal = mmsTimeTerm(tEndEff) + mmsReactionTerm(tEndEff);

        disp(' ')
        disp('###################################################')
        disp(['HOMOGENEOUS PRE-EQUILIBRATION. Step: ' num2str(solveCounter+1)])
        disp('###################################################')
        disp(' ')

        if ii == 1
            mmsForcingTermTotal_init = mmsTimeTerm(0) + mmsReactionTerm(0);
            homogeneousRxnSourceTerms_init = computeRHS_makrand_homogeneousPre(nSpecies, eqConcValues, constants);% + mmsForcingTermTotal_init;
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
            homogeneousRxnSourceTerms = computeRHS_makrand_homogeneousPre(nSpecies, eqConcValuesStar, constants);% + mmsForcingTermTotal;

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

            analyticalCharge = 0;
            for sIndex = 1:nSpecies
                analyticalCharge = analyticalCharge + constants.vale(1, sIndex) * mmsConc(t_rxn, sIndex)/constants.m3ToLiter;
            end

%             figure(111)
%             hold on
% %             plot(t_rxn, abs(eqConcValues(1) - mmsConc(t_rxn, 1)), 'ok')
% %             plot(t_rxn, abs(eqConcValues(2) - mmsConc(t_rxn, 2)), 'ob')
% %             plot(t_rxn, abs(eqConcValues(3) - mmsConc(t_rxn, 3)), 'or')
% %             plot(t_rxn, abs(eqConcValues(4) - mmsConc(t_rxn, 4)), 'om')
% %             plot(t_rxn, abs(eqConcValues(5) - mmsConc(t_rxn, 5)), 'og')
% %             plot(t_rxn, abs(eqConcValues(6) - mmsConc(t_rxn, 6)), 'oc')
% %             plot(t_rxn, abs(eqConcValues(7) - mmsConc(t_rxn, 7)), 'oy')
% 
% %             plot(t_rxn, eqConcValues(1) - mmsConc(t_rxn, 1), 'ok')
% %             plot(t_rxn, eqConcValues(2) - mmsConc(t_rxn, 2), 'ob')
% %             plot(t_rxn, eqConcValues(3) - mmsConc(t_rxn, 3), 'or')
% %             plot(t_rxn, eqConcValues(4) - mmsConc(t_rxn, 4), 'om')
% %             plot(t_rxn, eqConcValues(5) - mmsConc(t_rxn, 5), 'og')
% %             plot(t_rxn, eqConcValues(6) - mmsConc(t_rxn, 6), 'oc')
% %             plot(t_rxn, eqConcValues(7) - mmsConc(t_rxn, 7), 'oy')
% 
% %             plot(t_rxn, mmsConc(t_rxn, 1), 'ok')
% %             plot(t_rxn, mmsConc(t_rxn, 2), 'ob')
% %             plot(t_rxn, mmsConc(t_rxn, 3), 'or')
% %             plot(t_rxn, mmsConc(t_rxn, 4), 'om')
% %             plot(t_rxn, mmsConc(t_rxn, 5), 'og')
% %             plot(t_rxn, mmsConc(t_rxn, 6), 'oc')
% %             plot(t_rxn, mmsConc(t_rxn, 7), 'oy')
% % 
%             plot(t_rxn, eqConcValues(1), 'sk')
%             plot(t_rxn, eqConcValues(2), 'sb')
%             plot(t_rxn, eqConcValues(3), 'sr')
%             plot(t_rxn, eqConcValues(4), 'sm')
%             plot(t_rxn, eqConcValues(5), 'sg')
%             plot(t_rxn, eqConcValues(6), 'sc')
%             plot(t_rxn, eqConcValues(7), 'sy')
%             hold off
% 
%             figure(112)
%             hold on
%             plot(t_rxn, sum(constants.vale(1,:) .* eqConcValues')/constants.m3ToLiter, 'ob')
% % %             plot(t_rxn, analyticalCharge, 'sb')
%             hold off
%             xlabel('time (s)')
%             ylabel('net charge (M)')
        end
    end
end

% eqConcValues(1) = constants.initVal(1, 1) * constants.m3ToLiter;
constants.initVal = repmat(eqConcValues', [nX, 1]);
layerInfo(1).yRightBCValues = eqConcValues';

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

bulkElectrolyteC = zeros(nSpecies, 1);
for sIndex = 1:nSpecies
    bulkElectrolyteC(sIndex) = constants.initVal(1,sIndex) * constants.m3ToLiter;
end

% % Change concentration units to mol/m^3, to match diffusivity.
% for yIndex = 1:size(fields3D, 2)
%     for zIndex = 1:size(fields3D, 3)
%         for ii = 1:nSpecies
%             fields3D(ii:(nSpecies+1):end, yIndex, zIndex) = fields3D(ii:(nSpecies+1):end, yIndex, zIndex) * constants.m3ToLiter;
%         end
%     end
% end

% Create 3Dd-slice-form of BC vector for Dirichlet BC, in y direction
yLeftBCVec3DForm = zeros(size(fields3D(:, 1, :)));
for speciesIndex = 1:nSpecies
    layerInfo(1).yLeftBCValues(speciesIndex) = constants.m3ToLiter*layerInfo(1).yLeftBCValues(speciesIndex);
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
    layerInfo(1).yRightBCValues(speciesIndex) = constants.m3ToLiter*layerInfo(1).yRightBCValues(speciesIndex);
    for zIndex = 1:nZ
        if layerInfo(1).yRightBulkDirichlet
            yRightBCVec3DForm(speciesIndex:(nSpecies+1):end, 1, zIndex) = constants.initVal(1, speciesIndex);
        else
            yRightBCVec3DForm(speciesIndex:(nSpecies+1):end, 1, zIndex) = layerInfo(1).yRightBCValues(speciesIndex);
        end
    end
end
% Neumann BC for potential
yRightBCVec3DForm((nSpecies+1):(nSpecies+1):end, 1, :) = layerInfo(1).yRightPotentialBCValue;

% Create 3Dd-slice-form of BC vector for Dirichlet BC, in z direction
zBottomBCVec3DForm = zeros(size(fields3D(:, :, 1)));
for speciesIndex = 1:nSpecies
    layerInfo(1).zBottomBCValues(speciesIndex) = constants.m3ToLiter*layerInfo(1).zBottomBCValues(speciesIndex);
    for yIndex = 1:nY
        zBottomBCVec3DForm(speciesIndex:(nSpecies+1):end, yIndex, 1) = layerInfo(1).zBottomBCValues(speciesIndex);
    end
end
% Neumann BC for potential
zBottomBCVec3DForm((nSpecies+1):(nSpecies+1):end, :, 1) = layerInfo(1).zBottomPotentialBCValue;

zTopBCVec3DForm = zeros(size(fields3D(:, :, 1)));
for speciesIndex = 1:nSpecies
    layerInfo(1).zTopBCValues(speciesIndex) = constants.m3ToLiter*layerInfo(1).zTopBCValues(speciesIndex);
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
gasVelocity(1) = gasLayerInfo.gasFlowInlet / ...
    (gasLayerInfo.gasXYArea) / gasLayerInfo.gasTortuosity;

% Output frequency
outputInt = max(floor(dtOut/dt), 1);
plotInt = max(floor(dtPlot/dt), 1);

doPlot = plotCheck;
doOutputFunction = 1;       % 1 is on, will produce binary files
showUnconvergedPlots = 0;   % Show plots when time refinement is required

%output generation function here

if doOutputFunction
    doOutput(xCenter, yCenter, zCenter, uniqueSpecies, fields3D, gasFields, gasVelocity, gasLayerInfo.gasSpeciesNames, t, folderName);
end

% Plotting function here
if doPlot
    doAllPlots(fields3D((nSpecies+1)+1:end-(nSpecies+1), :, :), gasFields, gasVelocity, xCenter(2:end-1), yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, gasLayerInfo.gasSpeciesNames, constants, gasLayerInfo)
end

% Convenience variable initialization
totalCellsYZPlane = length(yCenter) * length(zCenter);
deltaFieldsForAccumulation = cell(totalCellsYZPlane, 1);
globalChargeComparison = 0;

% Initialize matrix for velocity size
VelLHSOperator = [[1 zeros(1, nZ-1); diag(1./dzC(1:end-1), 1) - diag(1./dzC, 0)], [zeros(nZ,1); 1/dzC(end)]];

% Save original exchange current density, for linear ramp (to avoid
% transient)
originalExchangeCurrents = zeros(length(faradaicRxnInfo), 1);
for rxnIndex = 1:length(faradaicRxnInfo)
    originalExchangeCurrents(rxnIndex) = faradaicRxnInfo(rxnIndex).exchangeCurrentDensity;
end

dvC= zeros(nX, nY, nZ);
dxC = xFace(2:end) - xFace(1:end-1);
dyC = yFace(2:end) - yFace(1:end-1);
dzC = zFace(2:end) - zFace(1:end-1);
for xIndex =1:nX
    for yIndex = 1:nY
        for zIndex = 1:nZ
            dvC(xIndex, yIndex, zIndex) = dxC(xIndex) * dyC(yIndex) * dzC(zIndex);
        end
    end
end

load('mmsFunctions.mat');
% MMS control flags
Period = 1e-3;
epsilon = 0.5;
coeff_time = 1;
coeff_diff = 1;
coeff_migr = 1;
coeff_hrxn = 1;
coeff_frxn = 1;

xField = repmat(reshape(xCenter(2:end-1), [nX-2,1]), [1, nY, nZ]);
yField = repmat(reshape(yCenter, [1,nY,1]), [nX-2, 1, nZ]);
zField = repmat(reshape(zCenter, [1,1,nZ]), [nX-2, nY, 1]);

c1Eq = eqConcValues(1);
c2Eq = eqConcValues(2);
c3Eq = eqConcValues(3);
c4Eq = eqConcValues(4);
c5Eq = eqConcValues(5);
c6Eq = eqConcValues(6);
c7Eq = eqConcValues(7);
phiEq = layerInfo(1).yRightPotentialBCValue;

D1 = constants.diff(1,1);
D2 = constants.diff(1,2);
D3 = constants.diff(1,3);
D4 = constants.diff(1,4);
D5 = constants.diff(1,5);
D6 = constants.diff(1,6);
D7 = constants.diff(1,7);

z1 = constants.vale(1,1);
z2 = constants.vale(1,2);
z3 = constants.vale(1,3);
z4 = constants.vale(1,4);
z5 = constants.vale(1,5);
z6 = constants.vale(1,6);
z7 = constants.vale(1,7);

poro = constants.poro(1);
tort = constants.tort(1);

T = constants.T;
kb = constants.kb;
na = constants.NA;
e = constants.e;

x = xField;
y = yField;
z = zField;
Lx = xFace(end) - xFace(1);
Ly = yFace(end) - yFace(1);
Lz = zFace(end) - zFace(1);

rxns = constants.rxns{1};
for rIndex = 1:size(rxns, 1)
    eval(strcat("hRxnRate", strtrim(num2str(rIndex)), " = rxns(rIndex, 5);"));
end

fRxn1ActivationEnergy = faradaicRxnInfo(1).activationEnergy;
fRxn1CRef = faradaicRxnInfo(1).cRef;
fRxn1ExchangeCurrent = faradaicRxnInfo(1).exchangeCurrentDensity;
fRxn1MoleProduct2 = faradaicRxnInfo(1).products(2);
fRxn1MoleReactant = faradaicRxnInfo(1).reactants;
fRxn1StandardElectrodePotential = faradaicRxnInfo(1).standardElectrodePotential;
fRxn1SurfaceRoughness = faradaicRxnInfo(1).surfaceRoughness;
fRxn1TransferCoefficient = faradaicRxnInfo(1).transferCoefficient;

fRxn2ActivationEnergy = faradaicRxnInfo(2).activationEnergy;
fRxn2ExchangeCurrent = faradaicRxnInfo(2).exchangeCurrentDensity;
fRxn2MoleProduct2 = faradaicRxnInfo(2).products(2);
fRxn2StandardElectrodePotential = faradaicRxnInfo(2).standardElectrodePotential;
fRxn2SurfaceRoughness = faradaicRxnInfo(2).surfaceRoughness;
fRxn2TransferCoefficient = faradaicRxnInfo(2).transferCoefficient;

phiElectrode = constants.phiElectrode;

mmsForcing = zeros(nX-2, nY, nZ, 8);
mmsFields = zeros(nX-2, nY, nZ, 8);
mmsNorm = zeros(totalSteps, nSpecies+2);

%% Time stepping loop starts here

solveCounter = 0;

% History term initialization
dtEffOld = dt;
dtEffOldOld = dtEffOld;

fields3DOld = fields3D;
fields3DOldOld = fields3DOld;
fields3DOldOldOld = fields3DOldOld;

yLeftBCVec3DFormOld = yLeftBCVec3DForm;

stepCtr = 1;

% constants.vale = zeros(nX, nSpecies);

% Time stepping loop
for ii = (ceil(t/dt)+1):totalSteps
    
    notConverged = 1;
    dtEff = dt;
%     if ii == 1
%         dtEff = 1e-6;
%     end

    targetTime = t+dt;
    
    % Modify reaction rates according to linear ramp (avoids nasty
    % transient that produces super thin zones near boundary, due to rapid
    % equilibration of chemical reactions)
    if targetTime <= layerInfo(1).reactionRampEndTime
        for rxnIndex = 1:length(faradaicRxnInfo)
            faradaicRxnInfo(rxnIndex).exchangeCurrentDensity = t/layerInfo(1).reactionRampEndTime * originalExchangeCurrents(rxnIndex);
        end
    end
    
    % Initial RHS terms are saved for Crank Nicolson startup
    if ii == 1
        [faradaicRxnSourceTerms_init, transportRHSTerms_init, homogeneousRxnSourceTerms_init, transverseRHSTerms_init] = computeRHS_makrand(dxC, dxF, dyC, dyF, dzC, dzF, ...
            nSpecies, fields3D, fields3DOld, layerInfo(1).xFrontNoFluxCondition, layerInfo(end).xBackNoFluxCondition, layerInfo(1).yLeftNoFluxCondition, layerInfo(1).yRightNoFluxCondition, layerInfo(1).zBottomNoFluxCondition, layerInfo(1).zTopNoFluxCondition, ...
            constants, yLeftBCVec3DFormOld, yRightBCVec3DForm, zBottomBCVec3DForm, zTopBCVec3DForm, faradaicRxnInfo, layerInfo(1).isNeutralGas, layerInfo(1).yRightElectrolyteOverride, bulkElectrolyteL, bulkElectrolyteC, layerInfo(1).potentialIterations, layerInfo(1).debugMode);
    end

    while (notConverged) || (t ~= targetTime)
        tic
        tEndEff = t + dtEff;
        
%         % Voltage sweep code is not functional for 3D, needs to be repaired if
%         % desired
%         if layerInfo(1).voltageSweepOn
%             if tEndEff <= timeVals(end)
%                 yRightBCVec3DForm((nSpecies+1):(nSpecies+1):end) = interp1(timeVals, voltageVals, tEndEff);
%                 layerInfo(1).yRightPotentialBCValue = interp1(timeVals, voltageVals, tEndEff);
%             else
%                 yRightBCVec3DForm((nSpecies+1):(nSpecies+1):end) = voltageVals(end);
%                 layerInfo(1).yRightPotentialBCValue = voltageVals(end);
%             end
%         end

        if stepCtr == 1
            time = (tEndEff - tStart)/2;
        else
            time = tEndEff;
        end

        mmsForcing(:, :, :, 1) = eval(strjoin(["force_c1_f(" strjoin(string(args1), ", ") ")"]));
        mmsForcing(:, :, :, 2) = eval(strjoin(["force_c2_f(" strjoin(string(args2), ", ") ")"]));
        mmsForcing(:, :, :, 3) = eval(strjoin(["force_c3_f(" strjoin(string(args3), ", ") ")"]));
        mmsForcing(:, :, :, 4) = eval(strjoin(["force_c4_f(" strjoin(string(args4), ", ") ")"]));
        mmsForcing(:, :, :, 5) = eval(strjoin(["force_c5_f(" strjoin(string(args5), ", ") ")"]));
        mmsForcing(:, :, :, 6) = eval(strjoin(["force_c6_f(" strjoin(string(args6), ", ") ")"]));
        mmsForcing(:, :, :, 7) = eval(strjoin(["force_c7_f(" strjoin(string(args7), ", ") ")"]));
        mmsForcing(:, :, :, 8) = eval(strjoin(["force_phi_f(" strjoin(string(argsphi), ", ") ")"]));
        
        mmsFields(:, :, :, 1) = eval(strjoin(["c1_f(" strjoin(string(cargs1), ", ") ")"]));
        mmsFields(:, :, :, 2) = eval(strjoin(["c2_f(" strjoin(string(cargs2), ", ") ")"]));
        mmsFields(:, :, :, 3) = eval(strjoin(["c3_f(" strjoin(string(cargs3), ", ") ")"]));
        mmsFields(:, :, :, 4) = eval(strjoin(["c4_f(" strjoin(string(cargs4), ", ") ")"]));
        mmsFields(:, :, :, 5) = eval(strjoin(["c5_f(" strjoin(string(cargs5), ", ") ")"]));
        mmsFields(:, :, :, 6) = eval(strjoin(["c6_f(" strjoin(string(cargs6), ", ") ")"]));
        mmsFields(:, :, :, 7) = eval(strjoin(["c7_f(" strjoin(string(cargs7), ", ") ")"]));
        mmsFields(:, :, :, 8) = eval(strjoin(["phi_f(" strjoin(string(cargsphi), ", ") ")"]));
        
        %%%%% SOLVE FOR GAS FIELD CONCENTRATIONS AND VELOCITY
        % Gas velocity ODE solve
        gasRHSCellCenters = zeros(nZ, 1);
        for sIndex = 1:nGasSpecies
            [speciesFaradaicReactionRate, speciesDissolutionRate] = computeGasProductionRate(gasFields, gasVelocity, fields3D, dxC, dxF, dyC, dyF, dzC, dzF, nX, nY, nZ, gasLayerInfo, constants, uniqueSpecies, faradaicRxnInfo, sIndex);
            gasRHSCellCenters = gasRHSCellCenters + (speciesFaradaicReactionRate + speciesDissolutionRate) * constants.NA * constants.kb * constants.T / (gasLayerInfo.gasPressure);
        end
        gasRHS = [gasVelocity(1); gasRHSCellCenters];
        
        gasVelocity = VelLHSOperator \ gasRHS;
        
        % RK4 method, explicit
        k1Frac = dtEff * rhsFunctionGasConc(gasFields, gasVelocity, fields3D, dxC, dxF, dyC, dyF, dzC, dzF, nX, nY, nZ, gasLayerInfo, constants, uniqueSpecies, faradaicRxnInfo);
        
        k2Frac = dtEff * rhsFunctionGasConc(gasFields + 1/2*k1Frac, gasVelocity, fields3D, dxC, dxF, dyC, dyF, dzC, dzF, nX, nY, nZ, gasLayerInfo, constants, uniqueSpecies, faradaicRxnInfo);
        
        k3Frac = dtEff * rhsFunctionGasConc(gasFields + 1/2*k2Frac, gasVelocity, fields3D, dxC, dxF, dyC, dyF, dzC, dzF, nX, nY, nZ, gasLayerInfo, constants, uniqueSpecies, faradaicRxnInfo);
        
        k4Frac = dtEff * rhsFunctionGasConc(gasFields + k3Frac, gasVelocity, fields3D, dxC, dxF, dyC, dyF, dzC, dzF, nX, nY, nZ, gasLayerInfo, constants, uniqueSpecies, faradaicRxnInfo);
        
        gasFields = gasFields ...
            + 1/6 * k1Frac ...
            + 1/3 * (k2Frac + k3Frac) ...
            + 1/6 * (k4Frac);
        
        % Henry's Law implementation for neutralGas species
        for sIndex = 1:nGasSpecies
            if layerInfo(1).isNeutralGas(sIndex)
                aqueousIndex = gasLayerInfo.isTrackedInElectrolyte(sIndex);
                for zIndex = 1:nZ
%                     fields3D(aqueousIndex, :, zIndex) = (gasLayerInfo.henrysLawConstant(sIndex) * gasFields((zIndex-1)*nGasSpecies + sIndex) * gasLayerInfo.gasPressure);
%                     fields3D(end-(nSpecies+1)+aqueousIndex, :, zIndex) = (gasLayerInfo.henrysLawConstant(sIndex) * gasFields((zIndex-1)*nGasSpecies + sIndex) * gasLayerInfo.gasPressure);
%                     yLeftBCVec3DForm(aqueousIndex:(nSpecies+1):end, 1, zIndex) = (gasLayerInfo.henrysLawConstant(sIndex) * gasFields((zIndex-1)*nGasSpecies + sIndex) * gasLayerInfo.gasPressure);
                    
%                     fields3D(aqueousIndex, :, zIndex) = min(tEndEff/0.5, 1e-5/0.5)*(gasLayerInfo.henrysLawConstant(sIndex) * gasFields((zIndex-1)*nGasSpecies + sIndex) * gasLayerInfo.gasPressure);
%                     fields3D(end-(nSpecies+1)+aqueousIndex, :, zIndex) = min(tEndEff/0.5, 1e-5/0.5)*(gasLayerInfo.henrysLawConstant(sIndex) * gasFields((zIndex-1)*nGasSpecies + sIndex) * gasLayerInfo.gasPressure);
%                     yLeftBCVec3DForm(aqueousIndex:(nSpecies+1):end, 1, zIndex) = min(tEndEff/0.5, 1e-5/0.5)*(gasLayerInfo.henrysLawConstant(sIndex) * gasFields((zIndex-1)*nGasSpecies + sIndex) * gasLayerInfo.gasPressure);

                    fields3D(aqueousIndex, :, zIndex) = eqConcValues(1);
                    fields3D(end-(nSpecies+1)+aqueousIndex, :, zIndex) = eqConcValues(1);
                    yLeftBCVec3DForm(aqueousIndex:(nSpecies+1):end, 1, zIndex) = eqConcValues(1);

                end
            end
        end

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Equation coupling, Faradaic Rxn. %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(' ')
        disp('######################################################')
        disp(['FARADAIC + TRANSPORT ITERATIONS. Step: ' num2str(solveCounter+1)])
        disp('######################################################')
        disp(' ')

        fields3DCoupled = fields3D;
        deltaFields = zeros(size(fields3D));
        
%         % Preparing variables for block iteration loop
%         couplingCharge = -inf;
%         oldCouplingCharge = inf;
%         couplingIteration = 0;
%         localChargeComparison = 0;
%         concNormInitial = 0;
%         isFirstStep = 1;
%         
%         % Block iteration to ensure Faradaic reactions do not produce
%         % spurious charge
%         while (abs(couplingCharge - oldCouplingCharge) > 1e-3*globalChargeComparison) && (abs(couplingCharge - oldCouplingCharge) > 1e-3*localChargeComparison) && couplingIteration < 3
%             couplingIteration = couplingIteration + 1;
% 
%             % Preparing variables for concentration solve
%             concNorm = inf;
%             oldConcNorm = -inf;
%             
%             concIterationCtr = 0;
%             isImproving = 1;
%             
%             % Concentration solve should converge in one step due to
%             % linearity, need to debug and see why this is not happening.
%             while concNorm > concNormInitial * 1e-6 && concIterationCtr < 10 && (abs(concNorm - oldConcNorm) > concNormInitial * 1e-9) && isImproving || concIterationCtr < 3
%                 % Update concentration fields without homogeneous reactions
%                 if concIterationCtr == 0
%                     [faradaicRxnSourceTerms, transportRHSTerms, ~, transverseTransportOriginal] = computeRHS_makrand(dxC, dxF, dyC, dyF, dzC, dzF, ...
%                         nSpecies, fields3DCoupled, fields3DOld, layerInfo(1).xFrontNoFluxCondition, layerInfo(end).xBackNoFluxCondition, layerInfo(1).yLeftNoFluxCondition, layerInfo(1).yRightNoFluxCondition, layerInfo(1).zBottomNoFluxCondition, layerInfo(1).zTopNoFluxCondition, ...
%                         constants, yLeftBCVec3DFormOld, yRightBCVec3DForm, zBottomBCVec3DForm, zTopBCVec3DForm, faradaicRxnInfo, layerInfo(1).isNeutralGas, layerInfo(1).yRightElectrolyteOverride, bulkElectrolyteL, bulkElectrolyteC, layerInfo(1).potentialIterations, layerInfo(1).debugMode);
%                 else
%                     [faradaicRxnSourceTerms, transportRHSTerms, ~, transverseTransportUpdated] = computeRHS_makrand(dxC, dxF, dyC, dyF, dzC, dzF, ...
%                         nSpecies, fields3DCoupled, fields3DOld, layerInfo(1).xFrontNoFluxCondition, layerInfo(end).xBackNoFluxCondition, layerInfo(1).yLeftNoFluxCondition, layerInfo(1).yRightNoFluxCondition, layerInfo(1).zBottomNoFluxCondition, layerInfo(1).zTopNoFluxCondition, ...
%                         constants, yLeftBCVec3DFormOld, yRightBCVec3DForm, zBottomBCVec3DForm, zTopBCVec3DForm, faradaicRxnInfo, layerInfo(1).isNeutralGas, layerInfo(1).yRightElectrolyteOverride, bulkElectrolyteL, bulkElectrolyteC, layerInfo(1).potentialIterations, layerInfo(1).debugMode);
%                 end
% 
% %                 if concIterationCtr > 0
% %                     disp(['Norm of change in transverse transport terms: ' num2str(norm(reshape(transverseTransportUpdated - transverseTransportOriginal, [nSpecies*(nX-2)*nY*nZ, 1])))])
% %                 end
%                 
%                 %%%%%% PARALLELIZE HERE FOR SPEED %%%%%%
%                 concResidual = zeros(nY, nZ);
%                 parfor indexYZPlane = 1:totalCellsYZPlane
%                     [yIndex, zIndex] = ind2sub([length(yCenter), length(zCenter)], indexYZPlane);
%                     
%                     xStripMiddle = fields3DCoupled(:, yIndex, zIndex);
%                     xStripOld = fields3DOld(:, yIndex, zIndex);
%                     xStripOldOld = fields3DOldOld(:, yIndex, zIndex);
%                     xStripOldOldOld = fields3DOldOldOld(:, yIndex, zIndex);
%                     
%                     if zIndex == 5 && yIndex == 5
%                         [deltaXStripOut, ~, residual, computedJacobian] = doTimeStep_concentrationCouplingMMS(dxC, dxF, nSpecies, xStripMiddle, xStripOld, xStripOldOld, xStripOldOldOld, layerInfo(1).xFrontNoFluxCondition, layerInfo(end).xBackNoFluxCondition, ...
%                             t, dtEff, dtEffOld, dtEffOldOld, constants, faradaicRxnInfo, layerInfo(1).isNeutralGas, faradaicRxnSourceTerms(:, yIndex, zIndex, :), faradaicRxnSourceTerms_init(:, yIndex, zIndex, :), transportRHSTerms(:, yIndex, zIndex, :), transportRHSTerms_init(:, yIndex, zIndex, :), stepCtr, layerInfo(1).debugMode, mmsForcing(:, yIndex, zIndex, :));
%                     else
%                         [deltaXStripOut, ~, residual, ~] = doTimeStep_concentrationCouplingMMS(dxC, dxF, nSpecies, xStripMiddle, xStripOld, xStripOldOld, xStripOldOldOld, layerInfo(1).xFrontNoFluxCondition, layerInfo(end).xBackNoFluxCondition, ...
%                             t, dtEff, dtEffOld, dtEffOldOld, constants, faradaicRxnInfo, layerInfo(1).isNeutralGas, faradaicRxnSourceTerms(:, yIndex, zIndex, :), faradaicRxnSourceTerms_init(:, yIndex, zIndex, :), transportRHSTerms(:, yIndex, zIndex, :), transportRHSTerms_init(:, yIndex, zIndex, :), stepCtr, layerInfo(1).debugMode, mmsForcing(:, yIndex, zIndex, :));
%                     end
%                     
%                     deltaFieldsForAccumulation{indexYZPlane} = deltaXStripOut;
%     
%                     concResidual(indexYZPlane) = sum(residual.^2);
% 
% %                     % Jacobian test
% %                     if yIndex == 5 && zIndex == 5 && concIterationCtr == 1 && couplingIteration == 1 && ii >= 100
% %                         [faradaicRHSTerms, transportRHSTerms, ~] = computeRHS_makrand(dxC, dxF, dyC, dyF, dzC, dzF, ...
% %                             nSpecies, fields3DCoupled, fields3DOld, layerInfo(1).xFrontNoFluxCondition, layerInfo(end).xBackNoFluxCondition, layerInfo(1).yLeftNoFluxCondition, layerInfo(1).yRightNoFluxCondition, layerInfo(1).zBottomNoFluxCondition, layerInfo(1).zTopNoFluxCondition, ...
% %                             constants, yLeftBCVec3DFormOld, yRightBCVec3DForm, zBottomBCVec3DForm, zTopBCVec3DForm, faradaicRxnInfo, layerInfo(1).isNeutralGas, layerInfo(1).yRightElectrolyteOverride, bulkElectrolyteL, bulkElectrolyteC, layerInfo(1).potentialIterations, layerInfo(1).debugMode);
% %                         bOriginal = zeros(nSpecies*(nX-2), 1);
% %                         for sIndex = 1:nSpecies
% %                             if stepCtr == 1
% %                                 bOriginal(sIndex:nSpecies:end) = bOriginal(sIndex:nSpecies:end) ...
% %                                     + 1/2*reshape(transportRHSTerms(:, 5, 5, sIndex) + faradaicRHSTerms(:, 5, 5, sIndex), [nX-2,1]) ...
% %                                     + 1/2*reshape(transportRHSTerms_init(:, 5, 5, sIndex) + faradaicRxnSourceTerms_init(:, 5, 5, sIndex), [nX-2,1]);
% % 
% %                                 bOriginal(sIndex:nSpecies:end) = bOriginal(sIndex:nSpecies:end) - constants.poro(1)/dtEff*(fields3DCoupled(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5) - fields3DOld(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5));
% % 
% %                             else
% %                                 bOriginal(sIndex:nSpecies:end) = bOriginal(sIndex:nSpecies:end) + reshape(transportRHSTerms(:, 5, 5, sIndex) + faradaicRHSTerms(:, 5, 5, sIndex), [nX-2,1]);
% % 
% %                                 timeGamma = 1/(dtEffOld * (1 + dtEffOld/dtEff));
% %                                 timeBeta = -timeGamma*(1 + dtEffOld/dtEff)^2;
% %                                 timeAlpha = -timeGamma - timeBeta;
% % 
% %                                 bOriginal(sIndex:nSpecies:end) = bOriginal(sIndex:nSpecies:end) - constants.poro(1)*(timeAlpha * fields3DCoupled(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5) ...
% %                                     + timeBeta * fields3DOld(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5) ...
% %                                     + timeGamma * fields3DOldOld(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5));
% %                             end
% %                         end
% % 
% %                         
% %                         
% %                         numericalJacobian = zeros((nX-2)*nSpecies, (nX-2)*nSpecies);
% %                         
% %                         for jxIndex = 2:nX-1
% %                             for jsIndex = 1:nSpecies
% % 
% %                                 epsilonVar = 1e-3 * eqConcValues(jsIndex);
% % 
% %                                 epsilon3D = zeros(nX*(nSpecies+1), nY, nZ);
% %                                 epsilon3D((jxIndex-1)*(nSpecies+1) + jsIndex, 5, 5) = epsilonVar;
% % 
% %                                 [faradaicRHSTerms, transportRHSTerms, ~] = computeRHS_makrand(dxC, dxF, dyC, dyF, dzC, dzF, ...
% %                                     nSpecies, fields3DCoupled + epsilon3D, fields3DOld, layerInfo(1).xFrontNoFluxCondition, layerInfo(end).xBackNoFluxCondition, layerInfo(1).yLeftNoFluxCondition, layerInfo(1).yRightNoFluxCondition, layerInfo(1).zBottomNoFluxCondition, layerInfo(1).zTopNoFluxCondition, ...
% %                                     constants, yLeftBCVec3DFormOld, yRightBCVec3DForm, zBottomBCVec3DForm, zTopBCVec3DForm, faradaicRxnInfo, layerInfo(1).isNeutralGas, layerInfo(1).yRightElectrolyteOverride, bulkElectrolyteL, bulkElectrolyteC, layerInfo(1).potentialIterations, layerInfo(1).debugMode);
% %                                 bDelta = zeros(nSpecies*(nX-2), 1);
% %                                 for sIndex = 1:nSpecies
% %                                     if stepCtr == 1
% %                                         bDelta(sIndex:nSpecies:end) = bDelta(sIndex:nSpecies:end) ...
% %                                             + 1/2*reshape(transportRHSTerms(:, 5, 5, sIndex) + faradaicRHSTerms(:, 5, 5, sIndex), [nX-2,1]) ...
% %                                             + 1/2*reshape(transportRHSTerms_init(:, 5, 5, sIndex) + faradaicRxnSourceTerms_init(:, 5, 5, sIndex), [nX-2,1]);
% % 
% %                                         bDelta(sIndex:nSpecies:end) = bDelta(sIndex:nSpecies:end) - constants.poro(1)/dtEff*((fields3DCoupled(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5) + epsilon3D(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5)) - fields3DOld(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5));
% % 
% %                                     else
% %                                         bDelta(sIndex:nSpecies:end) = bDelta(sIndex:nSpecies:end) + reshape(transportRHSTerms(:, 5, 5, sIndex) + faradaicRHSTerms(:, 5, 5, sIndex), [nX-2,1]);
% % 
% %                                         timeGamma = 1/(dtEffOld * (1 + dtEffOld/dtEff));
% %                                         timeBeta = -timeGamma*(1 + dtEffOld/dtEff)^2;
% %                                         timeAlpha = -timeGamma - timeBeta;
% % 
% %                                         bDelta(sIndex:nSpecies:end) = bDelta(sIndex:nSpecies:end) - constants.poro(1)*(timeAlpha * (fields3DCoupled(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5) + epsilon3D(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5)) ...
% %                                             + timeBeta * fields3DOld(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5) ...
% %                                             + timeGamma * fields3DOldOld(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5));
% % 
% %                                     end
% %                                 end
% % 
% %                                 numericalJacobian(:, (jxIndex-2)*(nSpecies) + jsIndex) = -(bDelta-bOriginal)/epsilonVar;
% %                             end
% %                         end
% % 
% %                         test = numericalJacobian - computedJacobian;
% %                         disp(['Norm of Jacobian difference: ' num2str(norm(numericalJacobian - computedJacobian)/norm(full(computedJacobian)))])
% %                     end
% 
%                 end
% 
% %                 if concIterationCtr > 0
% %                     disp(['Norm of change in compueted Jacobian: ' num2str(norm(full(computedJacobian - originalComputedJacobian)))])
% %                 end
%     
%                 oldConcNorm = concNorm;
%                 concNorm = sqrt(sum(concResidual, [1, 2]));
%                 
%                 if layerInfo(1).debugMode
%                     disp(['Faradaic Iteration ' num2str(concIterationCtr) ' conc. |residual|: ' num2str(concNorm, '%.4e')]);
%                 end
% 
%                 if isFirstStep
%                     concNormInitial = concNorm;
%                     isFirstStep = 0;
%                 end
% 
%                 % Update concentration fields
%                 for indexYZPlane = 1:totalCellsYZPlane
%                     [yIndex, zIndex] = ind2sub([length(yCenter), length(zCenter)], indexYZPlane);
%                     for sIndex = 1:nSpecies
%                         deltaFields(sIndex:nSpecies+1:end, yIndex, zIndex) = ...
%                             [0; deltaFieldsForAccumulation{indexYZPlane}(sIndex:nSpecies:end); 0];
%                     end
%                 end
%                 fields3DCoupled = fields3DCoupled + deltaFields;
% 
%                 if concNorm >= oldConcNorm
%                     isImproving = 0;
%                 end
% 
%                 concIterationCtr = concIterationCtr + 1;
% 
%                 for sIndex = 1:nSpecies
%                     if min(reshape(fields3DCoupled(sIndex:nSpecies+1:end, :, :), [nX*nY*nZ, 1])) < 0
%                         disp(['Minimum of ' uniqueSpecies{sIndex} ' is negative.'])
%                     end
%                 end
% 
%             end
%         
%             % Prepare variables for potential solve
%             potentialNorm = inf;
%             oldPotentialNorm = -inf;
%             potentialNormInitial = 0;
%             potentialIterationCtr = 0;
%             potentialMaxDelta = inf;
% 
% %             % Potential solve (nonlinear and 3D)
% %             while potentialNorm > potentialNormInitial * 1e-6 && potentialIterationCtr < 10 && (abs(potentialNorm - oldPotentialNorm) > potentialNormInitial * 1e-9) && potentialMaxDelta > eps(1)
% %                 
% %                 oldPotentialNorm = potentialNorm;
% %                 
% %                 [deltaPotential, potentialNorm] = doTimeStep_potentialCoupling(dxC, dxF, dyC, dyF, dzC, dzF, dtEff, dtEffOld, dtEffOldOld, ...
% %                     nSpecies, fields3DCoupled, fields3DOld, fields3DOldOld, fields3DOldOldOld ,layerInfo(1).xFrontNoFluxCondition, layerInfo(end).xBackNoFluxCondition, layerInfo(1).yLeftNoFluxCondition, layerInfo(1).yRightNoFluxCondition, layerInfo(1).zBottomNoFluxCondition, layerInfo(1).zTopNoFluxCondition, ...
% %                     constants, yLeftBCVec3DForm, yRightBCVec3DForm, zBottomBCVec3DForm, zTopBCVec3DForm, faradaicRxnInfo, layerInfo(1).isNeutralGas, layerInfo(1).yRightElectrolyteOverride, bulkElectrolyteL, bulkElectrolyteC, layerInfo(1).potentialIterations, layerInfo(1).debugMode, stepCtr);
% % 
% %                 fields3DCoupled(nSpecies+1:nSpecies+1:end, :, :) = ...
% %                     fields3DCoupled(nSpecies+1:nSpecies+1:end, :, :) + deltaPotential;
% % 
% %                 if potentialIterationCtr == 0
% %                     potentialNormInitial = potentialNorm;
% %                 end
% % 
% %                 potentialMaxDelta = max(abs(deltaPotential(:)));
% %                 if layerInfo(1).debugMode
% %                     disp(['Faradaic Iteration ' num2str(potentialIterationCtr) ' potential |residual|: ' num2str(potentialNorm, '%.4e')]);
% %                 end
% % 
% %                 potentialIterationCtr = potentialIterationCtr + 1;
% %             end
%             
%             % Compute net charge for checking purposes
%             netCharge = zeros(nX-2, nY, nZ);
%             for sIndex = 1:nSpecies
%                 netCharge = netCharge + constants.vale(1, sIndex) * (fields3DCoupled(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), :, :) - fields3D(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), :, :));
%             end
%             
%             oldCouplingCharge = couplingCharge;
%             [~, maxInd] = max(abs(netCharge(:)/constants.m3ToLiter));
%             couplingCharge = netCharge(maxInd) / constants.m3ToLiter;
%             localChargeComparison = abs(couplingCharge);
%             
%             integratedCharge = sum(netCharge .* dvC(2:end-1, :, :), [1 2 3]) / constants.m3ToLiter / (xFace(end)-xFace(1)) / (yFace(end)-yFace(1)) / (zFace(end)-zFace(1));
%             
%             if layerInfo(1).debugMode
%                 disp(' ')
%                 disp(['Coupling block iteration: ' num2str(couplingIteration) ', Net charge generated: ' num2str(integratedCharge)])
%                 disp(' ')
%             end
% 
%         end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%   Homogeneous and Faradaic Rxn.  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(' ')
        disp('######################################################')
        disp(['[HOMOGENEOUS EQUILIBRATION ITERATIONS. Step: ' num2str(solveCounter+1)])
        disp('#######################################################')
        disp(' ')

        isFirstStep = 1;
        concNormInitial = 0;
        for homogeneousIteration = 1:2
            % Preparing variables for concentration solve
            concNorm = inf;
            oldConcNorm = -inf;
            concNormInitial = 0;
            isFirstStep = 1;
            concIterationCtr = 0;
            isImproving = 1;
            
            % Concentration solve should converge in one step due to
            % linearity, need to debug and see why this is not happening.
            while (concNorm > concNormInitial * 1e-9 && concIterationCtr < 10 && (abs(concNorm - oldConcNorm) > concNormInitial * 1e-12) && isImproving) || concIterationCtr < 3
                % Update concentration fields without homogeneous reactions
                [faradaicRxnSourceTerms, transportRHSTerms, homogeneousRxnSourceTerms, ~] = computeRHS_makrand(dxC, dxF, dyC, dyF, dzC, dzF, ...
                    nSpecies, fields3DCoupled, fields3DOld, layerInfo(1).xFrontNoFluxCondition, layerInfo(end).xBackNoFluxCondition, layerInfo(1).yLeftNoFluxCondition, layerInfo(1).yRightNoFluxCondition, layerInfo(1).zBottomNoFluxCondition, layerInfo(1).zTopNoFluxCondition, ...
                    constants, yLeftBCVec3DFormOld, yRightBCVec3DForm, zBottomBCVec3DForm, zTopBCVec3DForm, faradaicRxnInfo, layerInfo(1).isNeutralGas, layerInfo(1).yRightElectrolyteOverride, bulkElectrolyteL, bulkElectrolyteC, layerInfo(1).potentialIterations, layerInfo(1).debugMode);
                
                %%%%%% PARALLELIZE HERE FOR SPEED %%%%%%
                concResidual = zeros(nY, nZ);
                for indexYZPlane = 1:totalCellsYZPlane
                    [yIndex, zIndex] = ind2sub([length(yCenter), length(zCenter)], indexYZPlane);
                    
                    xStripMiddle = fields3DCoupled(:, yIndex, zIndex);
                    xStripOld = fields3DOld(:, yIndex, zIndex);
                    xStripOldOld = fields3DOldOld(:, yIndex, zIndex);
                    xStripOldOldOld = fields3DOldOldOld(:, yIndex, zIndex);
                    
                    if yIndex == 5 && zIndex == 5
                        [deltaXStripOut, ~, residual, computedJacobian] = doTimeStep_concentrationWithHomogeneousMMS(dxC, dxF, nSpecies, xStripMiddle, xStripOld, xStripOldOld, xStripOldOldOld, layerInfo(1).xFrontNoFluxCondition, layerInfo(end).xBackNoFluxCondition, ...
                            t, dtEff, dtEffOld, dtEffOldOld, constants, faradaicRxnInfo, layerInfo(1).isNeutralGas, faradaicRxnSourceTerms(:, yIndex, zIndex, :), faradaicRxnSourceTerms_init(:, yIndex, zIndex, :), transportRHSTerms(:, yIndex, zIndex, :), transportRHSTerms_init(:, yIndex, zIndex, :), homogeneousRxnSourceTerms(:, yIndex, zIndex, :), homogeneousRxnSourceTerms_init(:, yIndex, zIndex, :), stepCtr, layerInfo(1).debugMode, mmsForcing(:, yIndex, zIndex, :));
                    else
                        [deltaXStripOut, ~, residual, ~] = doTimeStep_concentrationWithHomogeneousMMS(dxC, dxF, nSpecies, xStripMiddle, xStripOld, xStripOldOld, xStripOldOldOld, layerInfo(1).xFrontNoFluxCondition, layerInfo(end).xBackNoFluxCondition, ...
                            t, dtEff, dtEffOld, dtEffOldOld, constants, faradaicRxnInfo, layerInfo(1).isNeutralGas, faradaicRxnSourceTerms(:, yIndex, zIndex, :), faradaicRxnSourceTerms_init(:, yIndex, zIndex, :), transportRHSTerms(:, yIndex, zIndex, :), transportRHSTerms_init(:, yIndex, zIndex, :), homogeneousRxnSourceTerms(:, yIndex, zIndex, :), homogeneousRxnSourceTerms_init(:, yIndex, zIndex, :), stepCtr, layerInfo(1).debugMode, mmsForcing(:, yIndex, zIndex, :));
                    end
                    
                    deltaFieldsForAccumulation{indexYZPlane} = deltaXStripOut;
    
                    concResidual(indexYZPlane) = sum(residual.^2);


%                     % Jacobian test
%                     if yIndex == 5 && zIndex == 5 && homogeneousIteration > 1 && concIterationCtr == 1
%                         [faradaicRHSTerms, transportRHSTerms, homogeneousRHSTerms, ~] = computeRHS_makrand(dxC, dxF, dyC, dyF, dzC, dzF, ...
%                             nSpecies, fields3DCoupled, fields3DOld, layerInfo(1).xFrontNoFluxCondition, layerInfo(end).xBackNoFluxCondition, layerInfo(1).yLeftNoFluxCondition, layerInfo(1).yRightNoFluxCondition, layerInfo(1).zBottomNoFluxCondition, layerInfo(1).zTopNoFluxCondition, ...
%                             constants, yLeftBCVec3DFormOld, yRightBCVec3DForm, zBottomBCVec3DForm, zTopBCVec3DForm, faradaicRxnInfo, layerInfo(1).isNeutralGas, layerInfo(1).yRightElectrolyteOverride, bulkElectrolyteL, bulkElectrolyteC, layerInfo(1).potentialIterations, layerInfo(1).debugMode);
%                         bOriginal = zeros(nSpecies*(nX-2), 1);
%                         
%                         for sIndex = 1:nSpecies
%                             if stepCtr == 1
%                                 bOriginal(sIndex:nSpecies:end) = bOriginal(sIndex:nSpecies:end) ...
%                                     + 1/2*reshape(transportRHSTerms(:, 5, 5, sIndex) + faradaicRHSTerms(:, 5, 5, sIndex), [nX-2,1]) ...
%                                     + 1/2*reshape(transportRHSTerms_init(:, 5, 5, sIndex) + faradaicRxnSourceTerms_init(:, 5, 5, sIndex), [nX-2,1]) ...
%                                     + 1/2*reshape(homogeneousRxnSourceTerms_init(:, 5, 5, sIndex) + homogeneousRHSTerms(:, 5, 5, sIndex), [nX-2,1]);
% 
%                                 bOriginal(sIndex:nSpecies:end) = bOriginal(sIndex:nSpecies:end) - constants.poro(1)/dtEff*(fields3DCoupled(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5) - fields3DOld(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5));
% 
%                             else
%                                 bOriginal(sIndex:nSpecies:end) = bOriginal(sIndex:nSpecies:end) + reshape(transportRHSTerms(:, 5, 5, sIndex) + faradaicRHSTerms(:, 5, 5, sIndex) + homogeneousRHSTerms(:, 5, 5, sIndex), [nX-2,1]);
% 
%                                 timeGamma = 1/(dtEffOld * (1 + dtEffOld/dtEff));
%                                 timeBeta = -timeGamma*(1 + dtEffOld/dtEff)^2;
%                                 timeAlpha = -timeGamma - timeBeta;
% 
%                                 bOriginal(sIndex:nSpecies:end) = bOriginal(sIndex:nSpecies:end) - constants.poro(1)*(timeAlpha * fields3DCoupled(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5) ...
%                                     + timeBeta * fields3DOld(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5) ...
%                                     + timeGamma * fields3DOldOld(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5));
%                             end
%                         end
% 
% 
% 
%                         numericalJacobian = zeros((nX-2)*nSpecies, (nX-2)*nSpecies);
% 
%                         for jxIndex = 2:nX-1
%                             for jsIndex = 1:nSpecies
% 
%                                 epsilonVar = 1e-1 * eqConcValues(jsIndex);
% 
%                                 epsilon3D = zeros(nX*(nSpecies+1), nY, nZ);
%                                 epsilon3D((jxIndex-1)*(nSpecies+1) + jsIndex, 5, 5) = epsilonVar;
% 
%                                 [faradaicRHSTerms, transportRHSTerms, homogeneousRHSTerms, ~] = computeRHS_makrand(dxC, dxF, dyC, dyF, dzC, dzF, ...
%                                     nSpecies, fields3DCoupled + epsilon3D, fields3DOld, layerInfo(1).xFrontNoFluxCondition, layerInfo(end).xBackNoFluxCondition, layerInfo(1).yLeftNoFluxCondition, layerInfo(1).yRightNoFluxCondition, layerInfo(1).zBottomNoFluxCondition, layerInfo(1).zTopNoFluxCondition, ...
%                                     constants, yLeftBCVec3DFormOld, yRightBCVec3DForm, zBottomBCVec3DForm, zTopBCVec3DForm, faradaicRxnInfo, layerInfo(1).isNeutralGas, layerInfo(1).yRightElectrolyteOverride, bulkElectrolyteL, bulkElectrolyteC, layerInfo(1).potentialIterations, layerInfo(1).debugMode);
%                                 bDelta = zeros(nSpecies*(nX-2), 1);
%                                 for sIndex = 1:nSpecies
%                                     if stepCtr == 1
%                                         bDelta(sIndex:nSpecies:end) = bDelta(sIndex:nSpecies:end) ...
%                                             + 1/2*reshape(transportRHSTerms(:, 5, 5, sIndex) + faradaicRHSTerms(:, 5, 5, sIndex), [nX-2,1]) ...
%                                             + 1/2*reshape(transportRHSTerms_init(:, 5, 5, sIndex) + faradaicRxnSourceTerms_init(:, 5, 5, sIndex) - transverseRHSTerms_init(:, 5, 5, sIndex), [nX-2,1]) ...
%                                             + 1/2*reshape(homogeneousRxnSourceTerms_init(:, 5, 5, sIndex) + homogeneousRHSTerms(:, 5, 5, sIndex), [nX-2,1]);
% 
%                                         bDelta(sIndex:nSpecies:end) = bDelta(sIndex:nSpecies:end) - constants.poro(1)/dtEff*((fields3DCoupled(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5) + epsilon3D(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5)) - fields3DOld(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5));
% 
%                                     else
%                                         bDelta(sIndex:nSpecies:end) = bDelta(sIndex:nSpecies:end) + reshape(transportRHSTerms(:, 5, 5, sIndex) + faradaicRHSTerms(:, 5, 5, sIndex) + homogeneousRHSTerms(:, 5, 5, sIndex), [nX-2,1]);
% 
%                                         timeGamma = 1/(dtEffOld * (1 + dtEffOld/dtEff));
%                                         timeBeta = -timeGamma*(1 + dtEffOld/dtEff)^2;
%                                         timeAlpha = -timeGamma - timeBeta;
% 
%                                         bDelta(sIndex:nSpecies:end) = bDelta(sIndex:nSpecies:end) - constants.poro(1)*(timeAlpha * (fields3DCoupled(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5) + epsilon3D(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5)) ...
%                                             + timeBeta * fields3DOld(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5) ...
%                                             + timeGamma * fields3DOldOld(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), 5, 5));
% 
%                                     end
%                                 end
% 
%                                 numericalJacobian(:, (jxIndex-2)*(nSpecies) + jsIndex) = -(bDelta - bOriginal)/epsilonVar;
%                             end
%                         end
% 
%                         test = numericalJacobian - computedJacobian;
%                         disp(['Norm of Jacobian difference: ' num2str(norm(numericalJacobian - computedJacobian)/norm(full(computedJacobian)))])
%                     end


                end
                
                oldConcNorm = concNorm;
                concNorm = sqrt(sum(concResidual, [1, 2]));
                
                if layerInfo(1).debugMode
                    disp(['Homogeneous Iteration ' num2str(concIterationCtr) ' conc. |residual|: ' num2str(concNorm, '%.4e')]);
                end

                if isFirstStep
                    concNormInitial = concNorm;
                    isFirstStep = 0;
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
            while potentialNorm > potentialNormInitial * 1e-6 && potentialIterationCtr < 10 && (abs(potentialNorm - oldPotentialNorm) > potentialNormInitial * 1e-9) && potentialMaxDelta > eps(1)
                % Update electric potential
                oldPotentialNorm = potentialNorm;
                
                [deltaPotential, potentialNorm] = doTimeStep_potentialCouplingMMS(dxC, dxF, dyC, dyF, dzC, dzF, dtEff, dtEffOld, dtEffOldOld, ...
                    nSpecies, fields3DCoupled, fields3DOld, fields3DOldOld, fields3DOldOldOld ,layerInfo(1).xFrontNoFluxCondition, layerInfo(end).xBackNoFluxCondition, layerInfo(1).yLeftNoFluxCondition, layerInfo(1).yRightNoFluxCondition, layerInfo(1).zBottomNoFluxCondition, layerInfo(1).zTopNoFluxCondition, ...
                    constants, yLeftBCVec3DForm, yRightBCVec3DForm, zBottomBCVec3DForm, zTopBCVec3DForm, faradaicRxnInfo, layerInfo(1).isNeutralGas, layerInfo(1).yRightElectrolyteOverride, bulkElectrolyteL, bulkElectrolyteC, layerInfo(1).potentialIterations, layerInfo(1).debugMode, stepCtr, mmsForcing(:, :, :, 8));

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
                    netCharge = netCharge + constants.vale(1, sIndex) * (fields3DCoupled(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), :, :) - fields3D(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), :, :));
                end
                
                integratedCharge = sum(netCharge .* dvC(2:end-1, :, :), [1 2 3]) / constants.m3ToLiter / (xFace(end)-xFace(1)) / (yFace(end)-yFace(1)) / (zFace(end)-zFace(1));
                
                [~, maxInd] = max(abs(netCharge(:)/constants.m3ToLiter));
                homogeneousCharge = netCharge(maxInd) / constants.m3ToLiter;
                disp(' ')
                disp(['Homogeneos + Faradaic iteration ' num2str(homogeneousIteration) ', Net charge generated: ' num2str(integratedCharge)])
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
            netCharge = netCharge + constants.vale(1, sIndex) * (fields3DStar(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), :, :) - fields3D(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), :, :));
        end
        [maxDeltaChargeValue, maxInd] = max(abs(netCharge(:)/constants.m3ToLiter));
        
        bulkCharge = maxDeltaChargeValue > layerInfo(1).chargeThreshold;

        forceRefinement = 0;
        
        if nanVals || imagVals || forceRefinement% || bulkCharge || negVals
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
            dtEff = min([dtEff/2, 2*dtEffOld]);
            
            disp(['Refining time step to: ' num2str(dtEff)]);
            
            if doPlot && showUnconvergedPlots
                doAllPlots(fields3DStar, gasFields, gasVelocity, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, gasLayerInfo.gasSpeciesNames, constants, gasLayerInfo)
            end
            
        else
            t = tEndEff;

            notConverged = 0;
            solveCounter = solveCounter + 1;
            
            fields3D = fields3DStar;
            fields3DOldOldOld = fields3DOldOld;
            fields3DOldOld = fields3DOld;
            fields3DOld = fields3D;

            dtEffOldOld = dtEffOld;
            dtEffOld = dtEff;
            dtEff = min([targetTime - t, 1.1*dtEffOld]);

            yLeftBCVec3DFormOld = yLeftBCVec3DForm;

            netCharge = zeros(nX-2, nY, nZ);
            for sIndex = 1:nSpecies
                netCharge = netCharge + constants.vale(1, sIndex) * (fields3D(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), :, :));
            end
            [~, maxInd] = max(abs(netCharge(:)/constants.m3ToLiter));
            globalChargeComparison = abs(netCharge(maxInd)/constants.m3ToLiter);
            
            integratedCharge = sum(netCharge .* dvC(2:end-1, :, :), [1 2 3]) / constants.m3ToLiter / (xFace(end)-xFace(1)) / (yFace(end)-yFace(1)) / (zFace(end)-zFace(1));

            disp(' ')
            disp('---------------------------------------------------------')
            disp('---------------------------------------------------------')
            disp(['Reached time: ' num2str(t) ' using time step of ' num2str(dtEffOld) '.'])
            if layerInfo(1).debugMode
                disp(['Total charge at end of time step: ' num2str(integratedCharge)])
            end
            disp('---------------------------------------------------------')
            disp('---------------------------------------------------------')
            disp(' ')

            stepCtr = stepCtr + 1;

%             doAllPlots(fields3D, gasFields, gasVelocity, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, gasLayerInfo.gasSpeciesNames, constants, gasLayerInfo)
        end
        toc
    end
    
    t = targetTime;
    
    % Compute L2 error for MMS fields
    for sIndex = 1:nSpecies+1
        totalDifference = (mmsFields(:, :, :, sIndex) - fields3D(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), :, :)).^2;
        speciesNorm = (mmsFields(:, :, :, sIndex)).^2;
        mmsNorm(ii, sIndex+1) = sqrt(sum(totalDifference, [1 2 3])) / sqrt(sum(speciesNorm, [1 2 3]));
    end
    mmsNorm(ii, 1) = t;
    
    disp('')
    disp(['MMS norm ' uniqueSpecies{1} ': ' num2str(mmsNorm(ii, 2), '%.4e')])
    disp(['MMS norm ' uniqueSpecies{2} ': ' num2str(mmsNorm(ii, 3), '%.4e')])
    disp(['MMS norm ' uniqueSpecies{3} ': ' num2str(mmsNorm(ii, 4), '%.4e')])
    disp(['MMS norm ' uniqueSpecies{4} ': ' num2str(mmsNorm(ii, 5), '%.4e')])
    disp(['MMS norm ' uniqueSpecies{5} ': ' num2str(mmsNorm(ii, 6), '%.4e')])
    disp(['MMS norm ' uniqueSpecies{6} ': ' num2str(mmsNorm(ii, 7), '%.4e')])
    disp(['MMS norm ' uniqueSpecies{7} ': ' num2str(mmsNorm(ii, 8), '%.4e')])
    disp(['MMS norm potential: ' num2str(mmsNorm(ii, 9), '%.4e')])
    disp('')
    
     if (mod(ii, outputInt) == 0) && (doOutputFunction)
        % output generation function here
        doOutput(xCenter, yCenter, zCenter, uniqueSpecies, fields3D, gasFields, gasVelocity, gasLayerInfo.gasSpeciesNames, t, folderName);
        
        save([folderName 'mmsError.mat'], 'mmsNorm');
        
        disp(['Completed step ' num2str(ii) ' of ' num2str(totalSteps) ' (' num2str(100*ii/totalSteps) '%)'])
    end
    
    if (mod(ii, plotInt) == 0) && doPlot
        % Plotting function here
%         doAllPlots(fields3D, gasFields, gasVelocity, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, gasLayerInfo.gasSpeciesNames, constants, gasLayerInfo)
        
        figure(77)
        hold on
        plot(mmsNorm(1:ii, 1), mmsNorm(1:ii, 2), 'ok')
        plot(mmsNorm(1:ii, 1), mmsNorm(1:ii, 3), 'sb')
        plot(mmsNorm(1:ii, 1), mmsNorm(1:ii, 4), '^r')
        plot(mmsNorm(1:ii, 1), mmsNorm(1:ii, 5), 'vm')
        plot(mmsNorm(1:ii, 1), mmsNorm(1:ii, 6), 'dc')
        plot(mmsNorm(1:ii, 1), mmsNorm(1:ii, 7), 'xg')
        plot(mmsNorm(1:ii, 1), mmsNorm(1:ii, 8), 'or')
        hold off
        legend(uniqueSpecies)
        xlabel('time')
        set(gca, 'yscale', 'log')
        ylabel('||c - c_{mms}|| / ||c_{mms}||')
        drawnow
    end
end
end