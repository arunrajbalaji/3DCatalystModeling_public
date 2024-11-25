function [xStripMiddle, t] = doTimeStep_potentialFirst(dxC, dxF, dyC, dyFLeft, dyFRight, dzC, dzFBottom, dzFTop, ...
    nSpecies, xStripMiddle, xStripOld, xStripLeft, xStripRight, xStripBottom, xStripTop, ...
    isFrontFluxCondition, isBackFluxCondition, isLeftFluxCondition, isRightFluxCondition, isBottomFluxCondition, isTopFluxCondition, ...
    isFrontPotentialNeumannCondition, isBackPotentialNeumannCondition, isLeftPotentialNeumannCondition, isRightPotentialNeumannCondition, isBottomPotentialNeumannCondition, isTopPotentialNeumannCondition, ...
    t, dt, relaxationParameter, constants, rowInd, colInd, nnz, nEq, xFrontBCValues, xBackBCValues, yLeftBCValues, yRightBCValues, zBottomBCValues, zTopBCValues, ...
    xFrontPotentialBCValue, xBackPotentialBCValue, yLeftPotentialBCValue, yRightPotentialBCValue, zBottomPotentialBCValue, zTopPotentialBCValue, ...
    faradaicRxnInfo, isNeutralGas, yRightElectrolyteOverride, xStripLeftLeft, dyFLeftLeft, bulkElectrolyteL, bulkElectrolyteC, nIter)
% doTimeStep: Perform the time step using our first order implict scheme in
% time and second order central differences in space, with staggered mesh.
%
%   [z, t] = doTimeStep(xCenter, xFace, uniqueSpecies, z, t, dt, constants)
%
% Inputs:
%       dx_c            - Mesh size at cell centers.
%       dx_f            - Mesh size at cell faces.
%       nSpecies        - Number of uique species
%       z               - Unknowns vector before step is performed, ordered
%                         for speed (need to think about this)
%       t               - time before step is performed.
%       dt              - time step
%       constants       - physical constants, in struct as per output of
%                         genPhysConstArrays.m function.
%       rowInd          - LHS matrix A row indices, in correct order
%       colInd          - LHS matrix A column indices, in correct order
%
% Outputs:
%       z               - Unknows vector after times step
%       t               - Time after time step
%
% Example: 
%       Line 1 of example
%       Line 2 of example
%       Line 3 of example
%
% Other m-files required: none
% MAT-files required: none
%
% See also: genPhysConstArrays.m, makeLHSIndexVectors.m
%
% Author: Arunraj Balaji
% Stanford University, Mani Group
% email: abalaji@stanford.edu
% Last revision: 13-July-2021
%------------- BEGIN CODE --------------
% Unpack constants
poro = constants.poro;
tort = constants.tort;
perm = constants.perm;
diff = constants.diff;
acti = constants.acti;
vale = constants.vale;
e = constants.e;
kb = constants.kb;
T = constants.T;
NA = constants.NA;
m3ToLiter = constants.m3ToLiter;
backCharge = constants.backCharge;
wienBeta = constants.wienBeta;
phiElectrode = constants.phiElectrode;

effectiveDiffVec = diff(1,:) * poro(1) / tort(1);

maxWienFactor = 0;

% Might need this later, implemented here in a very ad-hoc way. Leaving for
% reference, definitely improve if necessary.
% deltaStern = 1e-10;
% exchCurrent = 0.0;
% convFactor = 1/(1e3 / 1e4 * 6.022e23 * 1.602e-19);
% faradaicCoeff = deltaStern*e/kb/T;

stericACubed = constants.stericACubed;
stericOnOffVec = constants.stericOnOffVec;

xN = length(dxC);

vals = zeros(nnz,1);
b = zeros(nEq,1);

t = t + dt;
dtEff = dt/nIter;

% yFluxRightDiffTotal = zeros(xN-2, 1);
% yFluxRightElectroTotal = zeros(xN-2, 1);
% yFluxRightStericTotal = zeros(xN-2, 1);
% yFluxRightAppliedTotal = zeros(xN-2, 1);
% 
% yFluxLeftDiffTotal = zeros(xN-2, 1);
% yFluxLeftElectroTotal = zeros(xN-2, 1);
% yFluxLeftStericTotal = zeros(xN-2, 1);
% yFluxLeftAppliedTotal = zeros(xN-2, 1);

xRHSTotal = zeros(xN-2, nSpecies);
yLeftDiffAllSpecies = zeros(xN-2, nSpecies);
yRightDiffAllSpecies = zeros(xN-2, nSpecies);
yLeftElectroAllSpecies = zeros(xN-2, nSpecies);
yRightElectroAllSpecies = zeros(xN-2, nSpecies);
yLeftStericAllSpecies = zeros(xN-2, nSpecies);
yRightStericAllSpecies = zeros(xN-2, nSpecies);
zBottomDiffAllSpecies = zeros(xN-2, nSpecies);
zTopDiffAllSpecies = zeros(xN-2, nSpecies);
zBottomElectroAllSpecies = zeros(xN-2, nSpecies);
zTopElectroAllSpecies = zeros(xN-2, nSpecies);
zBottomStericAllSpecies = zeros(xN-2, nSpecies);
zTopStericAllSpecies = zeros(xN-2, nSpecies);

% For use in electrolyte BC
xStripMiddleStar = xStripMiddle;

% disp('Start of middle iteration.')

% Loop over iterations
for iterIndex = 1:nIter
    
    xStripOld = xStripMiddle;
    
    % Must zero out b in advance (for reaction terms)
    b = zeros(nEq,1);
    
    for concentrationIteration = 1:1
    %% Determine LHS matrix A values and RHS vector b values
    ctr = 1;
    
    % Species equations
    % Loop over species
    for speciesIndex = 1:nSpecies
        
        % Loop over points, ignoring first and last cells
        % ColInd in matrix will be off by one set of variables compared to kk!!!
        for xLocIndex = 2:xN-1
            % Useful shortcuts, which only need to be computed once
            outerDerivativeFrontFace = 1/dxC(xLocIndex)*(poro(xLocIndex)+poro(xLocIndex-1))/2*min(diff(xLocIndex,speciesIndex),diff(xLocIndex-1,speciesIndex))/((tort(xLocIndex)+tort(xLocIndex-1))/2);
            outerDerivativeBackFace = 1/dxC(xLocIndex)*(poro(xLocIndex)+poro(xLocIndex+1))/2*min(diff(xLocIndex,speciesIndex),diff(xLocIndex+1,speciesIndex))/((tort(xLocIndex)+tort(xLocIndex+1))/2);
            gammaInvLeft = 2/(acti(xLocIndex,speciesIndex) + acti(xLocIndex-1,speciesIndex));
            gammaInvRight = 2/(acti(xLocIndex,speciesIndex) + acti(xLocIndex+1,speciesIndex));
            stericOnOff = stericOnOffVec(speciesIndex);
            xLeftStericSum = sum(xStripMiddle((xLocIndex-2)*(nSpecies+1)+1:(xLocIndex-1)*(nSpecies+1)-1).*stericOnOffVec);
            xCenterStericSum = sum(xStripMiddle((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1).*stericOnOffVec);
            xRightStericSum = sum(xStripMiddle((xLocIndex)*(nSpecies+1)+1:(xLocIndex+1)*(nSpecies+1)-1).*stericOnOffVec);
            
            % species cell to left
            if xLocIndex > 2
                leftElectromigrationTerm = nextLeftElectromigrationTerm;
                leftDiffTerm = nextLeftDiffTerm;
                leftStericTerm = nextLeftStericTerm;
                vals(ctr) = outerDerivativeFrontFace * (leftElectromigrationTerm - leftDiffTerm*acti(xLocIndex-1,speciesIndex) - leftStericTerm);
                ctr = ctr + 1;
            end
            
            % species cell to right
            if xLocIndex < xN-1
                rightElectromigrationTerm = 1/2*((vale(1,speciesIndex)*e/kb/T)*(xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1) + (nSpecies+1)) - xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1))))/dxF(xLocIndex);
                rightDiffTerm = 1/dxF(xLocIndex)*gammaInvRight;
                rightStericTerm = stericOnOff * 1/2 * (log(1 - stericACubed * xRightStericSum) ...
                - log(1 - stericACubed * xCenterStericSum))/dxF(xLocIndex);
                vals(ctr) = outerDerivativeBackFace * (-rightElectromigrationTerm - rightDiffTerm*acti(xLocIndex+1,speciesIndex) + rightStericTerm);
                nextLeftElectromigrationTerm = rightElectromigrationTerm;
                nextLeftDiffTerm = rightDiffTerm;
                nextLeftStericTerm = rightStericTerm;
                
                ctr = ctr + 1;
            end

            % species cell center
            if xLocIndex == 2
                if isFrontPotentialNeumannCondition
                    leftElectromigrationTerm = 1/2*(vale(1,speciesIndex)*e/kb/T)*xFrontPotentialBCValue;
                else
                    leftElectromigrationTerm = 1/2*((vale(1,speciesIndex)*e/kb/T)*(xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1)  + nSpecies + 1+(nSpecies+1)) - xStripMiddle(nSpecies+1)))/dxF(xLocIndex-1);
                end
                leftDiffTerm = 1/dxF(xLocIndex-1)*gammaInvLeft;
                leftStericTerm = stericOnOff * 1/2 * (log(1 - stericACubed * xCenterStericSum) ...
                - log(1 - stericACubed * xLeftStericSum))/dxF(xLocIndex-1);
            elseif xLocIndex == xN-1
                if isBackPotentialNeumannCondition
                    rightElectromigrationTerm = 1/2*(vale(1,speciesIndex)*e/kb/T)*xBackPotentialBCValue;
                else
                    rightElectromigrationTerm = 1/2*((vale(1,speciesIndex)*e/kb/T)*(xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1) + (nSpecies+1)) - xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1))))/dxF(xLocIndex);
                end
                rightDiffTerm = 1/dxF(xLocIndex)*gammaInvRight;
                rightStericTerm = stericOnOff * 1/2 * (log(1 - stericACubed * xRightStericSum) ...
                - log(1 - stericACubed * xCenterStericSum))/dxF(xLocIndex);
            end
            
            if xLocIndex == 2 && isFrontFluxCondition && ~isNeutralGas(speciesIndex)
                vals(ctr) = poro(xLocIndex)/dtEff ...
                    + outerDerivativeBackFace*(-rightElectromigrationTerm + rightDiffTerm*acti(xLocIndex,speciesIndex) + rightStericTerm);
            elseif xLocIndex == xN-1 && isBackFluxCondition && ~isNeutralGas(speciesIndex)
                vals(ctr) = poro(xLocIndex)/dtEff ...
                    + outerDerivativeFrontFace*(leftElectromigrationTerm + leftDiffTerm*acti(xLocIndex,speciesIndex) - leftStericTerm);
            else
                vals(ctr) = poro(xLocIndex)/dtEff ...
                    + outerDerivativeBackFace*(-rightElectromigrationTerm + rightDiffTerm*acti(xLocIndex,speciesIndex) + rightStericTerm)...
                    + outerDerivativeFrontFace*(leftElectromigrationTerm + leftDiffTerm*acti(xLocIndex,speciesIndex) - leftStericTerm);
            end
            ctr = ctr + 1;
            
            % Coupling to other species, for steric effects
            % Cell to left
            if xLocIndex > 2
                leftCrossSpecTerm = stericOnOff * (xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-2)*(nSpecies+1) + speciesIndex))/2 * ...
                    (stericACubed/(1 - stericACubed * xLeftStericSum))/dxF(xLocIndex-1);
                
                leftCrossSpecTermCenter = stericOnOff * (xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-2)*(nSpecies+1) + speciesIndex))/2 * ...
                    (stericACubed/(1 - stericACubed * xCenterStericSum))/dxF(xLocIndex-1);
                
                lStericCrossVal = -outerDerivativeFrontFace * leftCrossSpecTerm;
            end

            % Cell to right
            if xLocIndex < xN-1
                rightCrossSpecTerm = stericOnOff * (xStripMiddle((xLocIndex)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 * ...
                    (stericACubed/(1 - stericACubed * xRightStericSum))/dxF(xLocIndex);
                
                rightCrossSpecTermCenter = stericOnOff * (xStripMiddle((xLocIndex)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 * ...
                    (stericACubed/(1 - stericACubed * xCenterStericSum))/dxF(xLocIndex);
                
                rStericCrossVal = -outerDerivativeBackFace * rightCrossSpecTerm;
            end
            
            % Cell center
            if xLocIndex == 2
                leftCrossSpecTermCenter = stericOnOff * (xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-2)*(nSpecies+1) + speciesIndex))/2 * ...
                    (stericACubed/(1 - stericACubed * xCenterStericSum))/dxF(xLocIndex-1);
            elseif xLocIndex == xN-1
                rightCrossSpecTermCenter = stericOnOff * (xStripMiddle((xLocIndex)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 * ...
                    (stericACubed/(1 - stericACubed * xCenterStericSum))/dxF(xLocIndex);
            end
            centerStericCrossVal = outerDerivativeBackFace*rightCrossSpecTermCenter + outerDerivativeFrontFace*leftCrossSpecTermCenter;
            
            for rxnIndex=1:nSpecies
                % Cell to left
                if xLocIndex > 2
                    vals(ctr) = lStericCrossVal;
                    
                    ctr = ctr + 1;
                end
                
                % Cell to right
                if xLocIndex < xN-1
                    vals(ctr) = rStericCrossVal;
                    
                    ctr = ctr + 1;
                end
                
                % Cell center
                if xLocIndex == 2 && isFrontFluxCondition && ~isNeutralGas(speciesIndex)
                    vals(ctr) = outerDerivativeBackFace*rightCrossSpecTermCenter;
                elseif xLocIndex == xN-1 && isBackFluxCondition && ~isNeutralGas(speciesIndex)
                    vals(ctr) = outerDerivativeFrontFace*leftCrossSpecTermCenter;
                else
                    vals(ctr) = centerStericCrossVal;
                end
                ctr = ctr + 1;
            end
            
            if yRightElectrolyteOverride
                if (vale(1, speciesIndex) > 0) && (speciesIndex == 3)   % HYDROGEN
                    % Find OH- concentration at ghost cell center
                    ghostCellOH = -xStripLeft((xLocIndex-1)*(nSpecies+1) + 7) + 2*xStripMiddle((xLocIndex-1)*(nSpecies+1) + 7);
                    % Equilibrate water recombination
                    % Set Dirichlet value
                    K_eq = constants.rxns{1}(17,5)/constants.rxns{1}(18,5);
                    xStripRight((xLocIndex-1)*(nSpecies+1) + 7) = K_eq/ghostCellOH;
                elseif (vale(1, speciesIndex) > 0) && (speciesIndex == 6)   % POTASSIUM
                    % DO NOTHING (just use Dirichlet Eq. solution value)
                    
                    % Electroneutrality condition
%                     ghostCellCharge = sum(vale(1,:)'.* (-xStripLeft((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex-1)*(nSpecies+1)+nSpecies) + 2*xStripMiddle((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex-1)*(nSpecies+1)+nSpecies))) ...
%                         - vale(1,6)'.* (-xStripLeft((xLocIndex-1)*(nSpecies+1)+6) + 2*xStripMiddle((xLocIndex-1)*(nSpecies+1)+6))
%                     xStripRight((xLocIndex-1)*(nSpecies+1) + 6) = -ghostCellCharge;
                end
            end
            
            % RHS vector entries
            rightRHSDiffTerm = rightDiffTerm*(xStripMiddle((xLocIndex)*(nSpecies+1) + speciesIndex)*acti(xLocIndex+1,speciesIndex) - xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)*acti(xLocIndex,speciesIndex));
            leftRHSDiffTerm = leftDiffTerm*(xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)*acti(xLocIndex,speciesIndex) - xStripMiddle((xLocIndex-2)*(nSpecies+1) + speciesIndex)*acti(xLocIndex-1,speciesIndex));

            rightRHSStericTerm = -(xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex)*(nSpecies+1) + speciesIndex)) * rightStericTerm;
            leftRHSStericTerm = -(xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-2)*(nSpecies+1) + speciesIndex)) * leftStericTerm;
            
            effectiveDiff = diff(1, speciesIndex) * poro(xLocIndex) / tort(xLocIndex);

            %Y direction terms start here
            yFluxRightDiff = - effectiveDiff * (xStripRight((xLocIndex-1)*(nSpecies+1) + speciesIndex) - xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)) / dyFRight;
            yFluxLeftDiff = - effectiveDiff * (xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex) - xStripLeft((xLocIndex-1)*(nSpecies+1) + speciesIndex)) / dyFLeft;
            
            yFluxRightElectro = -vale(1,speciesIndex)*e/(kb*T)*effectiveDiff * (xStripRight((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 ...
                * (xStripRight((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)) - xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1))) / dyFRight;
            yFluxLeftElectro = -vale(1,speciesIndex)*e/(kb*T)*effectiveDiff * (xStripLeft((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 ...
                * (xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)) - xStripLeft((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1))) / dyFLeft;
            
            % Prepare for steric flux terms
            yStericSumRight = sum(xStripRight((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1).*stericOnOffVec);
            yStericSumLeft = sum(xStripLeft((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1).*stericOnOffVec);
            stericSumMiddle = sum(xStripMiddle((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1).*stericOnOffVec);
            
            yFluxRightSteric = -effectiveDiff * (xStripRight((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 ...
                * (log(1 - stericACubed * yStericSumRight) - log(1 - stericACubed * stericSumMiddle))/dyFRight;
            yFluxLeftSteric = -effectiveDiff * (xStripLeft((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 ...
                * (log(1 - stericACubed * stericSumMiddle) - log(1 - stericACubed * yStericSumLeft))/dyFLeft;
            
            % Z direction terms start here
            zFluxTopDiff = - effectiveDiff * (xStripTop((xLocIndex-1)*(nSpecies+1) + speciesIndex) - xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)) / dzFTop;
            zFluxBottomDiff = - effectiveDiff * (xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex) - xStripBottom((xLocIndex-1)*(nSpecies+1) + speciesIndex)) / dzFBottom;
            
            zFluxTopElectro = -vale(1,speciesIndex)*e/kb/T*effectiveDiff * (xStripTop((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 ...
                * (xStripTop((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)) - xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1))) / dzFTop;
            zFluxBottomElectro = -vale(1,speciesIndex)*e/kb/T*effectiveDiff * (xStripBottom((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 ...
                * (xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)) - xStripBottom((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1))) / dzFBottom;
            
            % Prepare for steric flux terms
            zStericSumTop = sum(xStripTop((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1).*stericOnOffVec);
            zStericSumBottom = sum(xStripBottom((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1).*stericOnOffVec);
            
            zFluxTopSteric = -effectiveDiff * (xStripTop((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 ...
                * (log(1 - stericACubed * zStericSumTop) - log(1 - stericACubed * stericSumMiddle))/dzFTop;
            zFluxBottomSteric = -effectiveDiff * (xStripBottom((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 ...
                * (log(1 - stericACubed * stericSumMiddle) - log(1 - stericACubed * zStericSumBottom))/dzFBottom;
            
            yFluxLeftApplied = 0;
            if isLeftFluxCondition && ~isNeutralGas(speciesIndex)
                yFluxLeftDiff = 0;
                yFluxLeftElectro = 0;
                yFluxLeftSteric = 0;
                yFluxLeftApplied = yLeftBCValues(speciesIndex);
            elseif isLeftFluxCondition && isNeutralGas(speciesIndex)
               yFluxLeftSteric = 0;
            end
            
            if yRightElectrolyteOverride
                if vale(1, speciesIndex) < 0
                    delta_1_2 = dyC/2;
                    delta_3_2 = dyC/2 + dyFLeft;
                    delta_5_2 = dyC/2 + dyFLeft + dyFLeftLeft;
                    beta = (-1/(delta_5_2 - delta_1_2)*(delta_1_2^2 - delta_5_2^2)) / (-(delta_5_2 - delta_3_2)/(delta_5_2 - delta_1_2) * (delta_1_2^2 - delta_5_2^2) + delta_3_2^2 - delta_5_2^2);
                    alpha = 1/(delta_5_2 - delta_1_2) - beta * (delta_5_2 - delta_3_2)/(delta_5_2 - delta_1_2);
                    gamma = -alpha - beta;
                    yFluxRightDiff = - effectiveDiff * (alpha * xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex) + beta * xStripLeft((xLocIndex-1)*(nSpecies+1) + speciesIndex) + gamma * xStripLeftLeft((xLocIndex-1)*(nSpecies+1) + speciesIndex));
                    yFluxRightElectro = -vale(1,speciesIndex)*e/(kb*T)*effectiveDiff * (-1/2 * xStripLeft((xLocIndex-1)*(nSpecies+1) + speciesIndex) + 3/2 * xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)) ...
                        * (gamma * xStripLeftLeft((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)) + beta * xStripLeft((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)) + alpha * xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)));
                    effectiveBCConc = -1/2 * xStripLeft + 3/2 * xStripMiddle;
                    yStericSumRight = sum(effectiveBCConc((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1).*stericOnOffVec);
                    yFluxRightSteric = -effectiveDiff * (-1/2 * xStripLeft((xLocIndex-1)*(nSpecies+1) + speciesIndex) + 3/2 * xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)) ...
                        * (log(1 - stericACubed * yStericSumRight) - log(1 - stericACubed * stericSumMiddle))/dyFRight;
                elseif vale(1, speciesIndex) == 0
                    yFluxRightDiff = -poro(1)/tort(1)*(diff(1,1)+diff(1,4))/2 / bulkElectrolyteL * (bulkElectrolyteC(speciesIndex) - (-1/2 * xStripLeft((xLocIndex-1)*(nSpecies+1) + speciesIndex) + 3/2 * xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)));
                    yFluxRightElectro = 0;
                    yFluxRightSteric = 0;
                    
                elseif vale(1, speciesIndex) > 0
                    effectiveBCConc = -1/2 * xStripLeft + 3/2 * xStripMiddle;
                    yStericSumRight = sum(effectiveBCConc((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1).*stericOnOffVec);
                    yFluxRightSteric = -effectiveDiff * (-1/2 * xStripLeft((xLocIndex-1)*(nSpecies+1) + speciesIndex) + 3/2 * xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)) ...
                        * (log(1 - stericACubed * yStericSumRight) - log(1 - stericACubed * stericSumMiddle))/dyFRight;
                end
            end
            
            yFluxRightApplied = 0;
            if (isRightFluxCondition && ~yRightElectrolyteOverride) %|| (isRightFluxCondition && speciesIndex == 6)
                yFluxRightDiff = 0;
                yFluxRightElectro = 0;
                yFluxRightSteric = 0;
                yFluxRightApplied = yRightBCValues(speciesIndex);
            end
            
            zFluxBottomApplied = 0;
            if isBottomFluxCondition
                zFluxBottomDiff = 0;
                zFluxBottomElectro = 0;
                zFluxBottomSteric = 0;
                zFluxBottomApplied = zBottomBCValues(speciesIndex);
            end
            
            zFluxTopApplied = 0;
            if isTopFluxCondition
                zFluxTopDiff = 0;
                zFluxTopElectro = 0;
                zFluxTopSteric = 0;
                zFluxTopApplied = zTopBCValues(speciesIndex);
            end
            
            if xLocIndex == 2 && isFrontFluxCondition && ~isNeutralGas(speciesIndex)
                b(rowInd(ctr-1)) = -poro(xLocIndex)*(xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)-xStripOld((xLocIndex-1)*(nSpecies+1) + speciesIndex))/dtEff...
                    + outerDerivativeBackFace*(rightElectromigrationTerm*(xStripMiddle((xLocIndex)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)) + rightRHSDiffTerm  + rightRHSStericTerm) ...
                    - 1/dyC * (yFluxRightDiff - yFluxLeftDiff) ...
                    - 1/dyC * (yFluxRightElectro - yFluxLeftElectro) ...
                    - 1/dyC * (yFluxRightSteric - yFluxLeftSteric) ...
                    - 1/dyC * (yFluxRightApplied - yFluxLeftApplied) ...
                    - 1/dzC * (zFluxTopDiff - zFluxBottomDiff) ...
                    - 1/dzC * (zFluxTopElectro - zFluxBottomElectro) ...
                    - 1/dzC * (zFluxTopSteric - zFluxBottomSteric) ...
                    - 1/dzC * (zFluxTopApplied - zFluxBottomApplied) ...
                    + 1/dxC(2) * xFrontBCValues(speciesIndex);
                
                xRHSTotal(xLocIndex-1, speciesIndex) = outerDerivativeBackFace*(rightElectromigrationTerm*(xStripMiddle((xLocIndex)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)) + rightRHSDiffTerm  + rightRHSStericTerm) ...
                    + 1/dxC(2) * xFrontBCValues(speciesIndex);
                
            elseif xLocIndex == xN-1 && isBackFluxCondition && ~isNeutralGas(speciesIndex)
                b(rowInd(ctr-1)) = -poro(xLocIndex)*(xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)-xStripOld((xLocIndex-1)*(nSpecies+1) + speciesIndex))/dtEff...
                    - outerDerivativeFrontFace*(leftElectromigrationTerm*(xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-2)*(nSpecies+1) + speciesIndex)) + leftRHSDiffTerm + leftRHSStericTerm)...
                    - 1/dyC * (yFluxRightDiff - yFluxLeftDiff) ...
                    - 1/dyC * (yFluxRightElectro - yFluxLeftElectro) ...
                    - 1/dyC * (yFluxRightSteric - yFluxLeftSteric) ...
                    - 1/dyC * (yFluxRightApplied - yFluxLeftApplied) ...
                    - 1/dzC * (zFluxTopDiff - zFluxBottomDiff) ...
                    - 1/dzC * (zFluxTopElectro - zFluxBottomElectro) ...
                    - 1/dzC * (zFluxTopSteric - zFluxBottomSteric) ...
                    - 1/dzC * (zFluxTopApplied - zFluxBottomApplied) ...
                    - 1/dxC(end-1) * xBackBCValues(speciesIndex);
                
                xRHSTotal(xLocIndex-1, speciesIndex) = - outerDerivativeFrontFace*(leftElectromigrationTerm*(xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-2)*(nSpecies+1) + speciesIndex)) + leftRHSDiffTerm + leftRHSStericTerm) ...
                    - 1/dxC(end-1) * xBackBCValues(speciesIndex);
            else
                b(rowInd(ctr-1)) = -poro(xLocIndex)*(xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)-xStripOld((xLocIndex-1)*(nSpecies+1) + speciesIndex))/dtEff...
                    + outerDerivativeBackFace*(rightElectromigrationTerm*(xStripMiddle((xLocIndex)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)) + rightRHSDiffTerm + rightRHSStericTerm)...
                    - outerDerivativeFrontFace*(leftElectromigrationTerm*(xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-2)*(nSpecies+1) + speciesIndex)) + leftRHSDiffTerm + leftRHSStericTerm)...
                    - 1/dyC * (yFluxRightDiff - yFluxLeftDiff) ...
                    - 1/dyC * (yFluxRightElectro - yFluxLeftElectro) ...
                    - 1/dyC * (yFluxRightSteric - yFluxLeftSteric) ...
                    - 1/dyC * (yFluxRightApplied - yFluxLeftApplied) ...
                    - 1/dzC * (zFluxTopDiff - zFluxBottomDiff) ...
                    - 1/dzC * (zFluxTopElectro - zFluxBottomElectro) ...
                    - 1/dzC * (zFluxTopSteric - zFluxBottomSteric) ...
                    - 1/dzC * (zFluxTopApplied - zFluxBottomApplied);
                
                xRHSTotal(xLocIndex-1, speciesIndex) = outerDerivativeBackFace*(rightElectromigrationTerm*(xStripMiddle((xLocIndex)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)) + rightRHSDiffTerm + rightRHSStericTerm) ...
                    - outerDerivativeFrontFace*(leftElectromigrationTerm*(xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-2)*(nSpecies+1) + speciesIndex)) + leftRHSDiffTerm + leftRHSStericTerm);
            end

            yLeftDiffAllSpecies(xLocIndex-1, speciesIndex) = yFluxLeftDiff;
            yRightDiffAllSpecies(xLocIndex-1, speciesIndex) = yFluxRightDiff;
            yLeftElectroAllSpecies(xLocIndex-1, speciesIndex) = yFluxLeftElectro;
            yRightElectroAllSpecies(xLocIndex-1, speciesIndex) = yFluxRightElectro;
            yLeftStericAllSpecies(xLocIndex-1, speciesIndex) = yFluxLeftSteric;
            yRightStericAllSpecies(xLocIndex-1, speciesIndex) = yFluxRightSteric;
            zBottomDiffAllSpecies(xLocIndex-1, speciesIndex) = zFluxBottomDiff;
            zTopDiffAllSpecies(xLocIndex-1, speciesIndex) = zFluxTopDiff;
            zBottomElectroAllSpecies(xLocIndex-1, speciesIndex) = zFluxBottomElectro;
            zTopElectroAllSpecies(xLocIndex-1, speciesIndex) = zFluxTopElectro;
            zBottomStericAllSpecies(xLocIndex-1, speciesIndex) = zFluxBottomSteric;
            zTopStericAllSpecies(xLocIndex-1, speciesIndex) = zFluxTopSteric;
        
        end
        
        if ~(isreal(vals) && ~any(isnan(vals)) && isreal(b) && ~any(isnan(b)))
            xStripMiddle = NaN(size(xStripMiddle));
            return
        end
        
    end
    
    % Add reaction terms as necessary
    % Loop over grid points
    for xLocIndex = 2:xN-1
        % Add reaction terms for each cell
        rxns = constants.rxns{xLocIndex};
        
        for rxnIndex = 1:size(rxns, 1)
            % Reactants side first
            if ((rxns(rxnIndex,1) ~= 0) && (rxns(rxnIndex,2) ~= 0))

                vals(ctr) = poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
                ctr = ctr + 1;
                
                vals(ctr) = poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
                ctr = ctr + 1;
                
                vals(ctr) = poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
                ctr = ctr + 1;

                vals(ctr) = poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
                ctr = ctr + 1;

                b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,1)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,1))...
                    - poro(xLocIndex)*rxns(rxnIndex,5)...
                    * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))...
                    * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));

                b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,2)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,2))...
                    - poro(xLocIndex)*rxns(rxnIndex,5)...
                    * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))...
                    * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));

                if rxns(rxnIndex,3) ~= 0
                    vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
                    ctr = ctr + 1;

                    vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
                    ctr = ctr + 1;

                    b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3))...
                        + poro(xLocIndex)*rxns(rxnIndex,5)...
                        * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))...
                        * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
                end
                
                if rxns(rxnIndex,4) ~= 0
                    vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
                    ctr = ctr + 1;
                    
                    vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
                    ctr = ctr + 1;
                    
                    b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4))...
                        + poro(xLocIndex)*rxns(rxnIndex,5)...
                        * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))...
                        * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
                end
                
            elseif ((rxns(rxnIndex,1) ~= 0) && (rxns(rxnIndex,2) == 0))
                vals(ctr) = poro(xLocIndex)*rxns(rxnIndex,5);
                ctr = ctr + 1;
                
                b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,1)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,1))...
                    - poro(xLocIndex)*rxns(rxnIndex,5) * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
                
                
                if rxns(rxnIndex,3) ~= 0
                    vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5);
                    ctr = ctr + 1;
                    
                    b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3))...
                        + poro(xLocIndex)*rxns(rxnIndex,5) * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
                    
                end
                
                if rxns(rxnIndex,4) ~= 0
                    vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5);
                    ctr = ctr + 1;
                    
                    b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4))...
                        + poro(xLocIndex)*rxns(rxnIndex,5) * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
                    
                end
                
            elseif ((rxns(rxnIndex,1) == 0) && (rxns(rxnIndex,2) ~= 0))
                vals(ctr) = poro(xLocIndex)*rxns(rxnIndex,5);
                ctr = ctr + 1;
                
                b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,2)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,2))...
                    - poro(xLocIndex)*rxns(rxnIndex,5) * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
                
                if rxns(rxnIndex,3) ~= 0
                    vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5);
                    ctr = ctr + 1;
                    
                    b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3))...
                        + poro(xLocIndex)*rxns(rxnIndex,5) * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
                end
                
                if rxns(rxnIndex,4) ~= 0
                    vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5);
                    ctr = ctr + 1;
                    
                    b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4))...
                        + poro(xLocIndex)*rxns(rxnIndex,5) * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
                end
                
            else
                % Add constant term to RHS
                if rxns(rxnIndex,3) ~= 0
                    
                    b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3))...
                        + poro(xLocIndex)*rxns(rxnIndex,5);
                end
                
                if rxns(rxnIndex,4) ~= 0
                    
                    b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4))...
                        + poro(xLocIndex)*rxns(rxnIndex,5);
                end
            end
        end
        
        % Loop over Faradaic reactions and add to right and left side as
        % necessary.
        for fRxnIndex = 1:length(faradaicRxnInfo)
            if faradaicRxnInfo(fRxnIndex).reactants ~= 0
                b((xLocIndex-2)*nSpecies + faradaicRxnInfo(fRxnIndex).reactants) =  b((xLocIndex-2)*nSpecies + faradaicRxnInfo(fRxnIndex).reactants) ...
                    - 1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleReactantsPerMoleElectron*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                    * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                    * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                    * xStripMiddle((xLocIndex-1)*(nSpecies+1) + faradaicRxnInfo(fRxnIndex).reactants)/faradaicRxnInfo(fRxnIndex).cRef;
                vals(ctr) =  1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleReactantsPerMoleElectron*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                    * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                    * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                    * 1/faradaicRxnInfo(fRxnIndex).cRef;
                ctr = ctr + 1;
            end
            
            if faradaicRxnInfo(fRxnIndex).products(1) ~= 0
                if faradaicRxnInfo(fRxnIndex).reactants ~= 0
                    b((xLocIndex-2)*nSpecies + faradaicRxnInfo(fRxnIndex).products(1)) = b((xLocIndex-2)*nSpecies + faradaicRxnInfo(fRxnIndex).products(1)) ...
                        +  1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(1)*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                        * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                        * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                        * xStripMiddle((xLocIndex-1)*(nSpecies+1) + faradaicRxnInfo(fRxnIndex).reactants)/faradaicRxnInfo(fRxnIndex).cRef;
                    vals(ctr) = - 1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(1)*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                        * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                        * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                        * 1/faradaicRxnInfo(fRxnIndex).cRef;
                    ctr = ctr + 1;
                else
                    b((xLocIndex-2)*nSpecies + faradaicRxnInfo(fRxnIndex).products(1)) = b((xLocIndex-2)*nSpecies + faradaicRxnInfo(fRxnIndex).products(1)) ...
                        + 1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(1)*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                        * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                        * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential));
                end
            end
            
            if faradaicRxnInfo(fRxnIndex).products(2) ~= 0
                if faradaicRxnInfo(fRxnIndex).reactants ~= 0
                    b((xLocIndex-2)*nSpecies + faradaicRxnInfo(fRxnIndex).products(2)) = b((xLocIndex-2)*nSpecies + faradaicRxnInfo(fRxnIndex).products(2)) ...
                        + 1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(2)*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                        * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                        * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                        * xStripMiddle((xLocIndex-1)*(nSpecies+1) + faradaicRxnInfo(fRxnIndex).reactants)/faradaicRxnInfo(fRxnIndex).cRef;
                    vals(ctr) = -1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(2)*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                        * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                        * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                        * 1/faradaicRxnInfo(fRxnIndex).cRef;
                    ctr = ctr + 1;
                else
                    b((xLocIndex-2)*nSpecies + faradaicRxnInfo(fRxnIndex).products(2)) = b((xLocIndex-2)*nSpecies + faradaicRxnInfo(fRxnIndex).products(2)) ...
                        + 1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(2)*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                        * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                        * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential));
                end
            end
        end
    end
    
    %% Solve for the delta and update
    spparms('spumoni', 0)
%     warning('off', 'MATLAB:nearlySingularMatrix')
    if isreal(vals) && isreal(b) && ~any(isnan(vals)) && ~any(isnan(b))
        A = sparse(rowInd, colInd, vals, nEq, nEq, length(vals));
        dxStripConcentrationOnly = A\b;
    else
        xStripMiddle = NaN(size(xStripMiddle));
        return
    end
    
    % Update the starred variable
    concentrationIndexingVector = [logical(zeros(nSpecies+1, 1)); ...
        repmat(logical([ones(nSpecies, 1); 0]), [xN-2, 1]); ...
        logical(zeros(nSpecies+1, 1))];
    xStripMiddle(concentrationIndexingVector) ...
        = xStripMiddle(concentrationIndexingVector) + dxStripConcentrationOnly;
    
    % Exit immediately if first iteration returns imaginary values
    if ~isreal(dxStripConcentrationOnly) || any(isnan(dxStripConcentrationOnly))
        return
    end
    end
    
    for potentialIteration = 1:1
    % Storage for diaogonals of linear solve for current conservation
    dPhiCenterDiag = zeros(xN-2, 1);
    dPhiFrontDiag = zeros(xN-3, 1);
    dPhiBackDiag = zeros(xN-3, 1);
    
    % Create right side for current conservation equation
    currConsb = zeros(xN-2,1);
    
    % Current conservation equation, loop over grid points
    for xLocIndex = 2:xN-1
        
        electroSumCenter = (vale(1, :).^2.*effectiveDiffVec)*xStripMiddle((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1);

        % x-direction current conservation boundary handling
        % Note: extrapolation for C fields is performed onto the ghost cell
        % center, and then interpolated back onto the boundary face. This
        % is the same as just extrapolating onto the boundary face.
        if xLocIndex == 2
            if isFrontFluxCondition
                if isFrontPotentialNeumannCondition
                    electroSumFront = (vale(1, :).^2.*effectiveDiffVec)*(3/2*xStripMiddle((nSpecies+1)+1:2*(nSpecies+1)-1) - 1/2*xStripMiddle(2*(nSpecies+1)+1:3*(nSpecies+1)-1));
                else
                    electroSumFront = (vale(1, :).^2.*effectiveDiffVec)*(3/2*xStripMiddle((nSpecies+1)+1:2*(nSpecies+1)-1) - 1/2*xStripMiddle(2*(nSpecies+1)+1:3*(nSpecies+1)-1));
                end
            else
                if isFrontPotentialNeumannCondition
                    electroSumFront = (vale(1, :).^2.*effectiveDiffVec)*xStripMiddle((xLocIndex-2)*(nSpecies+1)+1:(xLocIndex-1)*(nSpecies+1)-1);
                else
                    electroSumFront = (vale(1, :).^2.*effectiveDiffVec)*xStripMiddle((xLocIndex-2)*(nSpecies+1)+1:(xLocIndex-1)*(nSpecies+1)-1);
                end
            end
        else
            electroSumFront = (vale(1, :).^2.*effectiveDiffVec)*xStripMiddle((xLocIndex-2)*(nSpecies+1)+1:(xLocIndex-1)*(nSpecies+1)-1);
        end
        
        if xLocIndex == xN-1
            if isBackFluxCondition
                if isBackPotentialNeumannCondition
                    electroSumBack = (vale(1, :).^2.*effectiveDiffVec)*(3/2*xStripMiddle((xN-1)*(nSpecies+1)-nSpecies:(xN-1)*(nSpecies+1)-1) - 1/2*xStripMiddle((xN-2)*(nSpecies+1)-nSpecies:(xN-2)*(nSpecies+1)-1));
                else
                    electroSumBack = (vale(1, :).^2.*effectiveDiffVec)*(3/2*xStripMiddle((xN-1)*(nSpecies+1)-nSpecies:(xN-1)*(nSpecies+1)-1) - 1/2*xStripMiddle((xN-2)*(nSpecies+1)-nSpecies:(xN-2)*(nSpecies+1)-1));
                end
            else
                if isBackPotentialNeumannCondition
                    electroSumBack = (vale(1, :).^2.*effectiveDiffVec)*xStripMiddle((xLocIndex)*(nSpecies+1)+1:(xLocIndex+1)*(nSpecies+1)-1);
                else
                    electroSumBack = (vale(1, :).^2.*effectiveDiffVec)*xStripMiddle((xLocIndex)*(nSpecies+1)+1:(xLocIndex+1)*(nSpecies+1)-1);
                end
            end
        else
            electroSumBack = (vale(1, :).^2.*effectiveDiffVec)*xStripMiddle((xLocIndex)*(nSpecies+1)+1:(xLocIndex+1)*(nSpecies+1)-1);
        end

        % dPhi center location (BC handling already performed)
        dPhiCenterDiag(xLocIndex-1) = -e^2/(kb*T)/dxC(xLocIndex) ...
                * (-(electroSumCenter + electroSumBack)/2/dxF(xLocIndex) ...
                - (electroSumCenter + electroSumFront)/2/dxF(xLocIndex-1));

        % dPhi front location
        if xLocIndex > 2
            dPhiFrontDiag(xLocIndex-2) = -e^2/(kb*T)/dxC(xLocIndex) ...
                * ((electroSumCenter + electroSumFront)/2/dxF(xLocIndex-1));
        end
        
        % dPhi back location
        if xLocIndex < xN-1
            dPhiBackDiag(xLocIndex-1) = -e^2/(kb*T)/dxC(xLocIndex) ...
                * ((electroSumCenter + electroSumBack)/2/dxF(xLocIndex));
        end
        
        % Loop over homogenized faradaic reactions
        % Note that exchangeCurrentDensity is a signed quantity, since we take
        % the exchange current to be signed depending on anode vs cathode
        homogenizedFaradaicReactions = 0;
        for fRxnIndex = 1:length(faradaicRxnInfo)
            if faradaicRxnInfo(fRxnIndex).reactants ~= 0
                individualFaradaicReactionTerm = 1/NA*faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                    * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                    * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                    * xStripMiddle((xLocIndex-1)*(nSpecies+1) + faradaicRxnInfo(fRxnIndex).reactants)/faradaicRxnInfo(fRxnIndex).cRef;
                homogenizedFaradaicReactions = homogenizedFaradaicReactions + individualFaradaicReactionTerm;
                    
            else
                individualFaradaicReactionTerm = 1/NA*faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                    * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                    * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential));
                homogenizedFaradaicReactions = homogenizedFaradaicReactions + individualFaradaicReactionTerm;
            end
            
            dPhiCenterDiag(xLocIndex-1) = dPhiCenterDiag(xLocIndex-1) ...
                - individualFaradaicReactionTerm * faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T);
        end
        
        currConsb(xLocIndex-1) = -1/dyC * sum(e * vale(1,:) .* (yRightDiffAllSpecies(xLocIndex-1,:) - yLeftDiffAllSpecies(xLocIndex-1,:) + yRightElectroAllSpecies(xLocIndex-1,:) - yLeftElectroAllSpecies(xLocIndex-1,:) + yRightStericAllSpecies(xLocIndex-1,:) - yLeftStericAllSpecies(xLocIndex-1,:)), [2]) ...
            - 1/dzC * sum(e * vale(1,:) .* (zTopDiffAllSpecies(xLocIndex-1,:) - zBottomDiffAllSpecies(xLocIndex-1,:) + zTopElectroAllSpecies(xLocIndex-1,:) - zBottomElectroAllSpecies(xLocIndex-1,:) + zTopStericAllSpecies(xLocIndex-1,:) - zBottomStericAllSpecies(xLocIndex-1,:)), [2]) ...
            + sum(e * vale(1,:) .* xRHSTotal(xLocIndex-1, :), [2]) ...
            + homogenizedFaradaicReactions;

% Add reactions here once confident in current conservation solve
    end
    
    % Construct and solve linear system for current conservation
    currConsA = spdiags([[dPhiFrontDiag; 0], dPhiCenterDiag, [0; dPhiBackDiag]], [-1 0 1], xN-2, xN-2);
    
    if isreal(currConsA) && isreal(currConsb) && ~any(isnan(currConsA(:))) && ~any(isnan(currConsb))
        dxStripPhiOnly = currConsA\currConsb;
    else
        xStripMiddle = NaN(size(xStripMiddle));
        return
    end
    
    % Update the starred variable
    phiIndexingVector = [logical(zeros(nSpecies+1, 1)); ...
        repmat(logical([zeros(nSpecies, 1); 1]), [xN-2, 1]); ...
        logical(zeros(nSpecies+1, 1))];
    xStripMiddle(phiIndexingVector) ...
        = xStripMiddle(phiIndexingVector) + relaxationParameter * dxStripPhiOnly;
%     xStripMiddle(phiIndexingVector) ...
%         = zeros(size(xStripMiddle(phiIndexingVector)));

%     if isRightFluxCondition
%         disp(['Max dPhi value: ', num2str(max(dxStripPhiOnly))])
%     end
    
    % Exit immediately if first iteration returns imaginary values
    if ~isreal(dxStripPhiOnly) || any(isnan(dxStripPhiOnly))
        return
    end
    end
    
end

%------------- END OF CODE --------------
