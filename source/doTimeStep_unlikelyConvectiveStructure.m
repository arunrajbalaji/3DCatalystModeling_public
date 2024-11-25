function [xStripMiddle, t] = doTimeStep(dxC, dxF, dyC, dyFLeft, dyFRight, dzC, dzFBottom, dzFTop, ...
    nSpecies, xStripMiddle, xStripMiddleTransverse, xStripLeft, xStripRight, xStripBottom, xStripTop, ...
    isFrontFluxCondition, isBackFluxCondition, isFrontConvectiveCondition, isBackConvectiveCondition, isLeftFluxCondition, isRightFluxCondition, isBottomFluxCondition, isTopFluxCondition, ...
    isFrontPotentialNeumannCondition, isBackPotentialNeumannCondition, isLeftPotentialNeumannCondition, isRightPotentialNeumannCondition, isBottomPotentialNeumannCondition, isTopPotentialNeumannCondition, ...
    t, dt, constants, rowInd, colInd, nnz, nEq, xFrontBCValues, xBackBCValues, yLeftBCValues, yRightBCValues, zBottomBCValues, zTopBCValues, ...
    xFrontPotentialBCValue, xBackPotentialBCValue, yLeftPotentialBCValue, yRightPotentialBCValue, zBottomPotentialBCValue, zTopPotentialBCValue, ...
    faradaicRxnInfo)
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

% Loop over iterations
nIter = 3;
xStripOriginalTime = xStripMiddle;
for iterIndex = 1:nIter
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
            
            % BC handling (convective override, then flux. Dirichlet is
            % default
            if isFrontConvectiveCondition && xLocIndex == 2
                if vale(1,speciesIndex) > 0
                    
                elseif vale(1,speciesIndex) < 0
                    leftElectromigrationTerm = ;
                    leftDiffTerm = ;
                    leftStericTerm = ;
                    vals(ctr) = poro(xLocIndex)/dt ...
                    + outerDerivativeBackFace*(-rightElectromigrationTerm + rightDiffTerm*acti(xLocIndex,speciesIndex) + rightStericTerm)...
                    + outerDerivativeFrontFace*(leftElectromigrationTerm + leftDiffTerm*acti(xLocIndex,speciesIndex) - leftStericTerm);
                else
                    
                end
            else
                if xLocIndex == 2 && isFrontFluxCondition
                    vals(ctr) = poro(xLocIndex)/dt ...
                        + outerDerivativeBackFace*(-rightElectromigrationTerm + rightDiffTerm*acti(xLocIndex,speciesIndex) + rightStericTerm);
                end
            end
            
            if isBackConvectiveCondition
                % STUFF
            else
                if xLocIndex == xN-1 && isBackFluxCondition
                    vals(ctr) = poro(xLocIndex)/dt ...
                        + outerDerivativeFrontFace*(leftElectromigrationTerm + leftDiffTerm*acti(xLocIndex,speciesIndex) - leftStericTerm);
                end
            end
            
            if (xLocIndex ~= 2 && xLocIndex ~= xN-1) || (xLocIndex==2 && isFrontFluxCondition == 0) || (xLocIndex==xN-1 && isBackFluxCondition == 0)
                vals(ctr) = poro(xLocIndex)/dt ...
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
                if xLocIndex == 2 && isFrontFluxCondition
                    vals(ctr) = outerDerivativeBackFace*rightCrossSpecTermCenter;
                elseif xLocIndex == xN-1 && isBackFluxCondition
                    vals(ctr) = outerDerivativeFrontFace*leftCrossSpecTermCenter;
                else
                    vals(ctr) = centerStericCrossVal;
                end
                ctr = ctr + 1;
            end
            
            % RHS vector entries
            rightRHSDiffTerm = rightDiffTerm*(xStripMiddle((xLocIndex)*(nSpecies+1) + speciesIndex)*acti(xLocIndex+1,speciesIndex) - xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)*acti(xLocIndex,speciesIndex));
            leftRHSDiffTerm = leftDiffTerm*(xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)*acti(xLocIndex,speciesIndex) - xStripMiddle((xLocIndex-2)*(nSpecies+1) + speciesIndex)*acti(xLocIndex-1,speciesIndex));

            rightRHSStericTerm = -(xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex)*(nSpecies+1) + speciesIndex)) * rightStericTerm;
            leftRHSStericTerm = -(xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddle((xLocIndex-2)*(nSpecies+1) + speciesIndex)) * leftStericTerm;
            
            effectiveDiff = diff(1, speciesIndex) * poro(xLocIndex) / tort(xLocIndex);

            %Y direction terms start here
            yFluxRightDiff = - effectiveDiff * (xStripRight((xLocIndex-1)*(nSpecies+1) + speciesIndex) - xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1) + speciesIndex)) / dyFRight;
            yFluxLeftDiff = - effectiveDiff * (xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1) + speciesIndex) - xStripLeft((xLocIndex-1)*(nSpecies+1) + speciesIndex)) / dyFLeft;
            
            yFluxRightElectro = -vale(1,speciesIndex)*effectiveDiff * (xStripRight((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 ...
                * (xStripRight((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)) - xStripMiddleTransverse((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1))) / dyFRight;
            yFluxLeftElectro = -vale(1,speciesIndex)*effectiveDiff * (xStripLeft((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 ...
                * (xStripMiddleTransverse((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)) - xStripLeft((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1))) / dyFLeft;
            
            % Prepare for steric flux terms
            yStericSumRight = sum(xStripRight((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1).*stericOnOffVec);
            yStericSumLeft = sum(xStripLeft((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1).*stericOnOffVec);
            stericSumMiddle = sum(xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1).*stericOnOffVec);
            
            yFluxRightSteric = -effectiveDiff * (xStripRight((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 ...
                * (log(1 - stericACubed * yStericSumRight) - log(1 - stericACubed * stericSumMiddle))/dyFRight;
            yFluxLeftSteric = -effectiveDiff * (xStripLeft((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 ...
                * (log(1 - stericACubed * stericSumMiddle) - log(1 - stericACubed * yStericSumLeft))/dyFLeft;
            
            % Z direction terms start here
            zFluxTopDiff = - effectiveDiff * (xStripTop((xLocIndex-1)*(nSpecies+1) + speciesIndex) - xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1) + speciesIndex)) / dzFTop;
            zFluxBottomDiff = - effectiveDiff * (xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1) + speciesIndex) - xStripBottom((xLocIndex-1)*(nSpecies+1) + speciesIndex)) / dzFBottom;
            
            zFluxTopElectro = -vale(1,speciesIndex)*effectiveDiff * (xStripTop((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 ...
                * (xStripTop((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)) - xStripMiddleTransverse((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1))) / dzFTop;
            zFluxBottomElectro = -vale(1,speciesIndex)*effectiveDiff * (xStripBottom((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 ...
                * (xStripMiddleTransverse((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)) - xStripBottom((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1))) / dzFBottom;
            
            % Prepare for steric flux terms
            zStericSumTop = sum(xStripTop((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1).*stericOnOffVec);
            zStericSumBottom = sum(xStripBottom((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1).*stericOnOffVec);
            
            zFluxTopSteric = -effectiveDiff * (xStripTop((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 ...
                * (log(1 - stericACubed * zStericSumTop) - log(1 - stericACubed * stericSumMiddle))/dzFTop;
            zFluxBottomSteric = -effectiveDiff * (xStripBottom((xLocIndex-1)*(nSpecies+1) + speciesIndex) + xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1) + speciesIndex))/2 ...
                * (log(1 - stericACubed * stericSumMiddle) - log(1 - stericACubed * zStericSumBottom))/dzFBottom;
            
            yFluxLeftApplied = 0;
            if isLeftFluxCondition
                yFluxLeftDiff = 0;
                yFluxLeftElectro = 0;
                yFluxLeftSteric = 0;
                yFluxLeftApplied = yLeftBCValues(speciesIndex);
            end
            
            yFluxRightApplied = 0;
            if isRightFluxCondition
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
            
            if xLocIndex == 2 && isFrontFluxCondition
                b(rowInd(ctr-1)) = -poro(xLocIndex)*(xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)-xStripOriginalTime((xLocIndex-1)*(nSpecies+1) + speciesIndex))/dt...
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
            elseif xLocIndex == xN-1 && isBackFluxCondition
                b(rowInd(ctr-1)) = -poro(xLocIndex)*(xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)-xStripOriginalTime((xLocIndex-1)*(nSpecies+1) + speciesIndex))/dt...
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
            else
                b(rowInd(ctr-1)) = -poro(xLocIndex)*(xStripMiddle((xLocIndex-1)*(nSpecies+1) + speciesIndex)-xStripOriginalTime((xLocIndex-1)*(nSpecies+1) + speciesIndex))/dt...
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
            end
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
                dPhidXStar = ((xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1) + (nSpecies+1)) + xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)))/2 ...
                - (xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)) + xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1)))/2)/dxC(xLocIndex);
                wienExponential = exp(rxns(rxnIndex,6)*wienBeta*abs(dPhidXStar));
                vals(ctr) = poro(xLocIndex)*rxns(rxnIndex,5)*wienExponential;
                ctr = ctr + 1;
                
                if wienExponential > maxWienFactor && rxns(rxnIndex,6)~=0
                    maxWienFactor = wienExponential;
                end
                
                b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,1)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,1))...
                    - poro(xLocIndex)*rxns(rxnIndex,5) * wienExponential * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
                
                if rxns(rxnIndex,6) ~= 0
                    
                    if xLocIndex == 2
                        vals(ctr) = wienExponential*wienBeta*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))*sign(dPhidXStar)/(2*dxC(xLocIndex));
                        ctr = ctr + 1;
                        
                    elseif xLocIndex==xN-1
                        vals(ctr) = -wienExponential*wienBeta*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))*sign(dPhidXStar)/(2*dxC(xLocIndex));
                        ctr = ctr + 1;
                        
                    else
                        vals(ctr) = wienExponential*wienBeta*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))*sign(dPhidXStar)/(2*dxC(xLocIndex));
                        ctr = ctr + 1;
                
                        vals(ctr) = -wienExponential*wienBeta*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))*sign(dPhidXStar)/(2*dxC(xLocIndex));
                        ctr = ctr + 1;
                        
                    end
                end
                
                if rxns(rxnIndex,3) ~= 0
                    vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5) * wienExponential;
                    ctr = ctr + 1;
                    
                    b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3))...
                        + poro(xLocIndex)*rxns(rxnIndex,5) * wienExponential * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
                    
                    if rxns(rxnIndex,6) ~= 0
                        if xLocIndex == 2
                            vals(ctr) = -wienExponential*wienBeta*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))*sign(dPhidXStar)/(2*dxC(xLocIndex));
                            ctr = ctr + 1;
                            
                        elseif xLocIndex==xN-1
                            vals(ctr) = wienExponential*wienBeta*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))*sign(dPhidXStar)/(2*dxC(xLocIndex));
                            ctr = ctr + 1;
                            
                        else
                            vals(ctr) = -wienExponential*wienBeta*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))*sign(dPhidXStar)/(2*dxC(xLocIndex));
                            ctr = ctr + 1;
                            
                            vals(ctr) = wienExponential*wienBeta*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))*sign(dPhidXStar)/(2*dxC(xLocIndex));
                            ctr = ctr + 1;
                            
                        end
                    end
                end
                
                if rxns(rxnIndex,4) ~= 0
                    vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5) * wienExponential;
                    ctr = ctr + 1;
                    
                    b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4))...
                        + poro(xLocIndex)*rxns(rxnIndex,5) * wienExponential * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
                    
                    if rxns(rxnIndex,6) ~= 0
                        if xLocIndex == 2
                            vals(ctr) = -wienExponential*wienBeta*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))*sign(dPhidXStar)/(2*dxC(xLocIndex));
                            ctr = ctr + 1;
                            
                        elseif xLocIndex==xN-1
                            vals(ctr) = wienExponential*wienBeta*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))*sign(dPhidXStar)/(2*dxC(xLocIndex));
                            ctr = ctr + 1;
                            
                        else
                            vals(ctr) = -wienExponential*wienBeta*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))*sign(dPhidXStar)/(2*dxC(xLocIndex));
                            ctr = ctr + 1;
                            
                            vals(ctr) = wienExponential*wienBeta*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))*sign(dPhidXStar)/(2*dxC(xLocIndex));
                            ctr = ctr + 1;
                            
                        end
                    end
                end
                
            elseif ((rxns(rxnIndex,1) == 0) && (rxns(rxnIndex,2) ~= 0))
                wienExponential = exp(rxns(rxnIndex,6)*wienBeta*abs(((xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1) + (nSpecies+1)) + xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)))/2 ...
                - (xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)) + xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1)))/2)/dxC(xLocIndex)));
                vals(ctr) = poro(xLocIndex)*rxns(rxnIndex,5) * wienExponential;
                ctr = ctr + 1;
                
                b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,2)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,2))...
                    - poro(xLocIndex)*rxns(rxnIndex,5) * wienExponential * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
                
                if rxns(rxnIndex,3) ~= 0
                    vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5) * wienExponential;
                    ctr = ctr + 1;
                    
                    b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3))...
                        + poro(xLocIndex)*rxns(rxnIndex,5) * wienExponential * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
                end
                
                if rxns(rxnIndex,4) ~= 0
                    vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5) * wienExponential;
                    ctr = ctr + 1;
                    
                    b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4))...
                        + poro(xLocIndex)*rxns(rxnIndex,5) * wienExponential * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
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
    
%     xStripMiddle(concentrationIndexingVector) ...
%         = ones(size(xStripMiddle(concentrationIndexingVector)));
    
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
        
        electroSumCenter = (vale(1, :).^2.*diff(xLocIndex, :))*xStripMiddle((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1);

        % x-direction current conservation boundary handling
        % Note: extrapolation for C fields is performed onto the ghost cell
        % center, and then interpolated back onto the boundary face. This
        % is the same as just extrapolating onto the boundary face.
        if xLocIndex == 2
            if isFrontFluxCondition
                if isFrontPotentialNeumannCondition
                    dCdXFront = (-xFrontBCValues'./diff(xLocIndex, :))' - vale(1, :)'*e/kb/T*xFrontPotentialBCValue.*(3/2*xStripMiddle((nSpecies+1)+1:2*(nSpecies+1)-1) - 1/2*xStripMiddle(2*(nSpecies+1)+1:3*(nSpecies+1)-1));
                    dPhidXFront = xFrontPotentialBCValue;
                    diffSumFrontFace = (vale(1, :).*diff(xLocIndex, :))*dCdXFront;
                    electroSumFront = (vale(1, :).^2.*diff(xLocIndex, :))*(3/2*xStripMiddle((nSpecies+1)+1:2*(nSpecies+1)-1) - 1/2*xStripMiddle(2*(nSpecies+1)+1:3*(nSpecies+1)-1));
                else
                    dCdXFront = (-xFrontBCValues'./diff(xLocIndex, :))' - vale(1, :)'*e/kb/T*(xStripMiddle(2*(nSpecies+1)) - xStripMiddle(nSpecies+1))/dxF(1).*(3/2*xStripMiddle((nSpecies+1)+1:2*(nSpecies+1)-1) - 1/2*xStripMiddle(2*(nSpecies+1)+1:3*(nSpecies+1)-1));
                    dPhidXFront = (xStripMiddle((xLocIndex)*(nSpecies+1)) - xStripMiddle((xLocIndex-1)*(nSpecies+1)))/dxF(xLocIndex-1);
                    diffSumFrontFace = (vale(1, :).*diff(xLocIndex, :))*dCdXFront;
                    electroSumFront = (vale(1, :).^2.*diff(xLocIndex, :))*(3/2*xStripMiddle((nSpecies+1)+1:2*(nSpecies+1)-1) - 1/2*xStripMiddle(2*(nSpecies+1)+1:3*(nSpecies+1)-1));
                end
            else
                if isFrontPotentialNeumannCondition
                    diffSumFrontFace = (vale(1, :).*diff(xLocIndex, :))*(xStripMiddle((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1)-xStripMiddle((xLocIndex-2)*(nSpecies+1)+1:(xLocIndex-1)*(nSpecies+1)-1))/dxF(xLocIndex-1);
                    dPhidXFront = xFrontPotentialBCValue;
                    electroSumFront = (vale(1, :).^2.*diff(xLocIndex, :))*xStripMiddle((xLocIndex-2)*(nSpecies+1)+1:(xLocIndex-1)*(nSpecies+1)-1);
                else
                    diffSumFrontFace = (vale(1, :).*diff(xLocIndex, :))*(xStripMiddle((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1)-xStripMiddle((xLocIndex-2)*(nSpecies+1)+1:(xLocIndex-1)*(nSpecies+1)-1))/dxF(xLocIndex-1);
                    dPhidXFront = (xStripMiddle((xLocIndex)*(nSpecies+1)) - xStripMiddle((xLocIndex-1)*(nSpecies+1)))/dxF(xLocIndex-1);
                    electroSumFront = (vale(1, :).^2.*diff(xLocIndex, :))*xStripMiddle((xLocIndex-2)*(nSpecies+1)+1:(xLocIndex-1)*(nSpecies+1)-1);
                end
            end
        else
            diffSumFrontFace = (vale(1, :).*diff(xLocIndex, :))*(xStripMiddle((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1)-xStripMiddle((xLocIndex-2)*(nSpecies+1)+1:(xLocIndex-1)*(nSpecies+1)-1))/dxF(xLocIndex-1);
            dPhidXFront = (xStripMiddle((xLocIndex)*(nSpecies+1)) - xStripMiddle((xLocIndex-1)*(nSpecies+1)))/dxF(xLocIndex-1);
            electroSumFront = (vale(1, :).^2.*diff(xLocIndex, :))*xStripMiddle((xLocIndex-2)*(nSpecies+1)+1:(xLocIndex-1)*(nSpecies+1)-1);
        end
        
        if xLocIndex == xN-1
            if isBackFluxCondition
                if isBackPotentialNeumannCondition
                    dCdXBack = (-xBackBCValues'./diff(xLocIndex, :))' - vale(1, :)'*e/kb/T*xBackPotentialBCValue.*(3/2*xStripMiddle((xN-1)*(nSpecies+1)-nSpecies:(xN-1)*(nSpecies+1)-1) - 1/2*xStripMiddle((xN-2)*(nSpecies+1)-nSpecies:(xN-2)*(nSpecies+1)-1));
                    dPhidXBack = xBackPotentialBCValue;
                    diffSumBackFace = (vale(1, :).*diff(xLocIndex, :))*dCdXBack;
                    electroSumBack = (vale(1, :).^2.*diff(xLocIndex, :))*(3/2*xStripMiddle((xN-1)*(nSpecies+1)-nSpecies:(xN-1)*(nSpecies+1)-1) - 1/2*xStripMiddle((xN-2)*(nSpecies+1)-nSpecies:(xN-2)*(nSpecies+1)-1));
                else
                    dCdXBack = (-xBackBCValues'./diff(xLocIndex, :))' - vale(1, :)'*e/kb/T*(xStripMiddle(end) - xStripMiddle(end-(nSpecies+1)))/dxF(end).*(3/2*xStripMiddle((xN-1)*(nSpecies+1)-nSpecies:(xN-1)*(nSpecies+1)-1) - 1/2*xStripMiddle((xN-2)*(nSpecies+1)-nSpecies:(xN-2)*(nSpecies+1)-1));
                    dPhidXBack = (xStripMiddle((xLocIndex+1)*(nSpecies+1)) - xStripMiddle((xLocIndex)*(nSpecies+1)))/dxF(xLocIndex);
                    diffSumBackFace = (vale(1, :).*diff(xLocIndex, :))*dCdXBack;
                    electroSumBack = (vale(1, :).^2.*diff(xLocIndex, :))*(3/2*xStripMiddle((xN-1)*(nSpecies+1)-nSpecies:(xN-1)*(nSpecies+1)-1) - 1/2*xStripMiddle((xN-2)*(nSpecies+1)-nSpecies:(xN-2)*(nSpecies+1)-1));
                end
            else
                if isBackPotentialNeumannCondition
                    diffSumBackFace = (vale(1, :).*diff(xLocIndex, :))*(xStripMiddle((xLocIndex)*(nSpecies+1)+1:(xLocIndex+1)*(nSpecies+1)-1)-xStripMiddle((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1))/dxF(xLocIndex);
                    dPhidXBack = xBackPotentialBCValue;
                    electroSumBack = (vale(1, :).^2.*diff(xLocIndex, :))*xStripMiddle((xLocIndex)*(nSpecies+1)+1:(xLocIndex+1)*(nSpecies+1)-1);
                else
                    diffSumBackFace = (vale(1, :).*diff(xLocIndex, :))*(xStripMiddle((xLocIndex)*(nSpecies+1)+1:(xLocIndex+1)*(nSpecies+1)-1)-xStripMiddle((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1))/dxF(xLocIndex);
                    dPhidXBack = (xStripMiddle((xLocIndex+1)*(nSpecies+1)) - xStripMiddle((xLocIndex)*(nSpecies+1)))/dxF(xLocIndex);
                    electroSumBack = (vale(1, :).^2.*diff(xLocIndex, :))*xStripMiddle((xLocIndex)*(nSpecies+1)+1:(xLocIndex+1)*(nSpecies+1)-1);
                end
            end
        else
            diffSumBackFace = (vale(1, :).*diff(xLocIndex, :))*(xStripMiddle((xLocIndex)*(nSpecies+1)+1:(xLocIndex+1)*(nSpecies+1)-1)-xStripMiddle((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1))/dxF(xLocIndex);
            dPhidXBack = (xStripMiddle((xLocIndex+1)*(nSpecies+1)) - xStripMiddle((xLocIndex)*(nSpecies+1)))/dxF(xLocIndex);
            electroSumBack = (vale(1, :).^2.*diff(xLocIndex, :))*xStripMiddle((xLocIndex)*(nSpecies+1)+1:(xLocIndex+1)*(nSpecies+1)-1);
        end
        
        % y-direction current conservation boundary handling
        % Note: extrapolation for C fields is performed onto the face, this
        % is the same is extrapolating onto an external center and then
        % interpolating back onto the boundary face.
        if isLeftFluxCondition
            if isLeftPotentialNeumannCondition
                dCdYLeft = (-yLeftBCValues./diff(xLocIndex, :))' - vale(1, :)'*e/kb/T*yLeftPotentialBCValue.*(3/2*xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:xLocIndex*(nSpecies+1)-1) - 1/2*xStripRight((xLocIndex-1)*(nSpecies+1)+1:xLocIndex*(nSpecies+1)-1));
                diffSumLeftFace = (vale(1, :).*diff(xLocIndex, :))*dCdYLeft;
                
                dPhidYLeft = yLeftPotentialBCValue;
                electroSumLeft = (vale(1, :).^2.*diff(xLocIndex, :))*(2*xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1) - xStripRight((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1));
            else
                dCdYLeft = (-yLeftBCValues./diff(xLocIndex, :))' - vale(1, :)'*e/kb/T*(xStripMiddleTransverse((xLocIndex)*(nSpecies+1)) - xStripLeft((xLocIndex)*(nSpecies+1)))/dyFLeft.*(3/2*xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:xLocIndex*(nSpecies+1)-1) - 1/2*xStripRight((xLocIndex-1)*(nSpecies+1)+1:xLocIndex*(nSpecies+1)-1));
                diffSumLeftFace = (vale(1, :).*diff(xLocIndex, :))*dCdYLeft;
                
                dPhidYLeft = (xStripMiddleTransverse((xLocIndex)*(nSpecies+1)) - xStripLeft((xLocIndex)*(nSpecies+1)))/dyFLeft;
                electroSumLeft = (vale(1, :).^2.*diff(xLocIndex, :))*(2*xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1) - xStripRight((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1));
            end
        else
            if isLeftPotentialNeumannCondition
                diffSumLeftFace = (vale(1, :).*diff(xLocIndex, :))*(xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1)-xStripLeft((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1))/dyFLeft;
                dPhidYLeft = yLeftPotentialBCValue;
                electroSumLeft = (vale(1, :).^2.*diff(xLocIndex, :))*xStripLeft((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1);
            else
                diffSumLeftFace = (vale(1, :).*diff(xLocIndex, :))*(xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1)-xStripLeft((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1))/dyFLeft;
                dPhidYLeft = (xStripMiddleTransverse((xLocIndex)*(nSpecies+1)) - xStripLeft((xLocIndex)*(nSpecies+1)))/dyFLeft;
                electroSumLeft = (vale(1, :).^2.*diff(xLocIndex, :))*xStripLeft((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1);
            end
        end
        
        if isRightFluxCondition
            if isRightPotentialNeumannCondition
                dCdYRight = (-yRightBCValues./diff(xLocIndex, :))' - vale(1, :)'*e/kb/T*yRightPotentialBCValue.*(3/2*xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:xLocIndex*(nSpecies+1)-1) - 1/2*xStripLeft((xLocIndex-1)*(nSpecies+1)+1:xLocIndex*(nSpecies+1)-1));
                diffSumRightFace = (vale(1, :).*diff(xLocIndex, :))*dCdYRight;
                
                dPhidYRight = yRightPotentialBCValue;
                electroSumRight = (vale(1, :).^2.*diff(xLocIndex, :))*(2*xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1) - xStripLeft((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1));
            else
                dCdYRight = (-yRightBCValues./diff(xLocIndex, :))' - vale(1, :)'*e/kb/T*(xStripRight((xLocIndex)*(nSpecies+1)) - xStripMiddleTransverse((xLocIndex)*(nSpecies+1)))/dyFRight.*(3/2*xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:xLocIndex*(nSpecies+1)-1) - 1/2*xStripLeft((xLocIndex-1)*(nSpecies+1)+1:xLocIndex*(nSpecies+1)-1));
                diffSumRightFace = (vale(1, :).*diff(xLocIndex, :))*dCdYRight;
                
                dPhidYRight = (xStripRight((xLocIndex)*(nSpecies+1)) - xStripMiddleTransverse((xLocIndex)*(nSpecies+1)))/dyFRight;
                electroSumRight = (vale(1, :).^2.*diff(xLocIndex, :))*(2*xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1) - xStripLeft((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1));
            end
        else
            if isRightPotentialNeumannCondition
                diffSumRightFace = (vale(1, :).*diff(xLocIndex, :))*(xStripRight((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1)-xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1))/dyFRight;
                dPhidYRight = yRightPotentialBCValue;
                electroSumRight = (vale(1, :).^2.*diff(xLocIndex, :))*xStripRight((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1);
            else
                diffSumRightFace = (vale(1, :).*diff(xLocIndex, :))*(xStripRight((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1)-xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1))/dyFRight;
                dPhidYRight = (xStripRight((xLocIndex)*(nSpecies+1)) - xStripMiddleTransverse((xLocIndex)*(nSpecies+1)))/dyFRight;
                electroSumRight = (vale(1, :).^2.*diff(xLocIndex, :))*xStripRight((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1);
            end
        end
        
        % z-direction current conservation boundary handling
        % Note: extrapolation for C fields is performed onto the face, this
        % is the same is extrapolating onto an external center and then
        % interpolating back onto the boundary face.
        if isBottomFluxCondition
            if isBottomPotentialNeumannCondition
                dCdZBottom = (-zBottomBCValues./diff(xLocIndex, :))' - vale(1, :)'*e/kb/T*zBottomPotentialBCValue.*(3/2*xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:xLocIndex*(nSpecies+1)-1) - 1/2*xStripTop((xLocIndex-1)*(nSpecies+1)+1:xLocIndex*(nSpecies+1)-1));
                diffSumBottomFace = (vale(1, :).*diff(xLocIndex, :))*dCdZBottom;
                
                dPhidZBottom = zBottomPotentialBCValue;
                electroSumBottom = (vale(1, :).^2.*diff(xLocIndex, :))*(2*xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1) - xStripTop((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1));
            else
                dCdZBottom = (-zBottomBCValues./diff(xLocIndex, :))' - vale(1, :)'*e/kb/T*(xStripMiddleTransverse((xLocIndex)*(nSpecies+1)) - xStripBottom((xLocIndex)*(nSpecies+1)))/dzFBottom.*(3/2*xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:xLocIndex*(nSpecies+1)-1) - 1/2*xStripTop((xLocIndex-1)*(nSpecies+1)+1:xLocIndex*(nSpecies+1)-1));
                diffSumBottomFace = (vale(1, :).*diff(xLocIndex, :))*dCdZBottom;
                
                dPhidZBottom = (xStripMiddleTransverse((xLocIndex)*(nSpecies+1)) - xStripBottom((xLocIndex)*(nSpecies+1)))/dzFBottom;
                electroSumBottom = (vale(1, :).^2.*diff(xLocIndex, :))*(2*xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1) - xStripTop((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1));
            end
        else
            if isBottomPotentialNeumannCondition
                diffSumBottomFace = (vale(1, :).*diff(xLocIndex, :))*(xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1)-xStripBottom((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1))/dzFBottom;
                dPhidZBottom = zBottomPotentialBCValue;
                electroSumBottom = (vale(1, :).^2.*diff(xLocIndex, :))*xStripBottom((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1);
            else
                diffSumBottomFace = (vale(1, :).*diff(xLocIndex, :))*(xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1)-xStripBottom((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1))/dzFBottom;
                dPhidZBottom = (xStripMiddleTransverse((xLocIndex)*(nSpecies+1)) - xStripBottom((xLocIndex)*(nSpecies+1)))/dzFBottom;
                electroSumBottom = (vale(1, :).^2.*diff(xLocIndex, :))*xStripBottom((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1);
            end
        end
        
        if isTopFluxCondition
            if isTopPotentialNeumannCondition
                dCdZTop = (-zTopBCValues./diff(xLocIndex, :))' - vale(1, :)'*e/kb/T*zTopPotentialBCValue.*(3/2*xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:xLocIndex*(nSpecies+1)-1) - 1/2*xStripBottom((xLocIndex-1)*(nSpecies+1)+1:xLocIndex*(nSpecies+1)-1));
                diffSumTopFace = (vale(1, :).*diff(xLocIndex, :))*dCdZTop;
                
                dPhidZTop = zTopPotentialBCValue;
                electroSumTop = (vale(1, :).^2.*diff(xLocIndex, :))*(2*xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1) - xStripBottom((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1));
            else
                dCdZTop = (-zTopBCValues./diff(xLocIndex, :))' - vale(1, :)'*e/kb/T*(xStripTop((xLocIndex)*(nSpecies+1)) - xStripMiddleTransverse((xLocIndex)*(nSpecies+1)))/dzFTop.*(3/2*xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:xLocIndex*(nSpecies+1)-1) - 1/2*xStripBottom((xLocIndex-1)*(nSpecies+1)+1:xLocIndex*(nSpecies+1)-1));
                diffSumTopFace = (vale(1, :).*diff(xLocIndex, :))*dCdZTop;
                
                dPhidZTop = (xStripTop((xLocIndex)*(nSpecies+1)) - xStripMiddleTransverse((xLocIndex)*(nSpecies+1)))/dzFTop;
                electroSumTop = (vale(1, :).^2.*diff(xLocIndex, :))*(2*xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1) - xStripBottom((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1));
            end
        else
            if isTopPotentialNeumannCondition
                diffSumTopFace = (vale(1, :).*diff(xLocIndex, :))*(xStripTop((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1)-xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1))/dzFTop;
                dPhidZTop = zTopPotentialBCValue;
                electroSumTop = (vale(1, :).^2.*diff(xLocIndex, :))*xStripTop((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1);
            else
                diffSumTopFace = (vale(1, :).*diff(xLocIndex, :))*(xStripTop((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1)-xStripMiddleTransverse((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1))/dzFTop;
                dPhidZTop = (xStripTop((xLocIndex)*(nSpecies+1)) - xStripMiddleTransverse((xLocIndex)*(nSpecies+1)))/dzFTop;
                electroSumTop = (vale(1, :).^2.*diff(xLocIndex, :))*xStripTop((xLocIndex-1)*(nSpecies+1)+1:(xLocIndex)*(nSpecies+1)-1);
            end
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
            homogenizedFaradaicReactions = homogenizedFaradaicReactions + ...
                1/NA*faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                * xStripMiddle((xLocIndex-1)*(nSpecies+1) + faradaicRxnInfo(fRxnIndex).reactants)/faradaicRxnInfo(fRxnIndex).cRef;
        end

        currConsb(xLocIndex-1) = - e/dxC(xLocIndex)*(diffSumBackFace - diffSumFrontFace) ...
            - e/dyC*(diffSumRightFace - diffSumLeftFace) ...            ...
            - e/dzC*(diffSumTopFace - diffSumBottomFace) ...
            + e^2/(kb*T)/dyC*((electroSumCenter + electroSumRight)/2*dPhidYRight ...
            - (electroSumCenter + electroSumLeft)/2*dPhidYLeft) ...
            + e^2/(kb*T)/dzC*((electroSumCenter + electroSumTop)/2*dPhidZTop ...
            - (electroSumCenter + electroSumBottom)/2*dPhidZBottom) ...
            + e^2/(kb*T)/dxC(xLocIndex-1)*((electroSumCenter + electroSumBack)/2*dPhidXBack ...
            - (electroSumCenter + electroSumFront)/2*dPhidXFront) ...
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
        = xStripMiddle(phiIndexingVector) + dxStripPhiOnly;
%     xStripMiddle(phiIndexingVector) ...
%         = zeros(size(xStripMiddle(phiIndexingVector)));
    
    % Exit immediately if first iteration returns imaginary values
    if ~isreal(dxStripPhiOnly) || any(isnan(dxStripPhiOnly))
        return
    end
    end
    
    
    %     for kk = 2:xN-1
%         % Poisson's equation
%         % potential cell to left
%         if kk > 2
%             leftPotentialTerm = nextLeftPotentialTerm;
%             vals(ctr) = leftPotentialTerm/dxC(kk);
%             ctr = ctr + 1;
%         end
%         
%         % potential cell to right
%         if kk < xN-1
%             rightPotentialTerm = 1/dxF(kk)*(perm(kk)+perm(kk+1))/2;
%             vals(ctr) = rightPotentialTerm/dxC(kk);
%             nextLeftPotentialTerm = rightPotentialTerm;
%             ctr = ctr + 1;
%         end
%         
%         % potential cell center
%         if kk == xN-1
%             rightPotentialTerm = 1/dxF(kk)*(perm(kk)+perm(kk+1))/2;
%         elseif kk == 2
%             leftPotentialTerm = 1/dxF(kk-1)*(perm(kk)+perm(kk-1))/2;
%         end
%         vals(ctr) = -leftPotentialTerm/dxC(kk) - rightPotentialTerm/dxC(kk);
%         ctr = ctr + 1;
%         
%         % Species cell centers
%         for jj = 1:nSpecies
%             vals(ctr) = NA*m3ToLiter*e*poro(kk)*vale(1,jj);
%             ctr = ctr + 1;
%             
%             if jj == 1
%                 % Reset RHS vector to zero before adding the new entries
%                 b(rowInd(ctr-1)) = -NA*m3ToLiter*poro(kk)*e*vale(1,jj)*xStripMiddle(rowInd(ctr-1)+jj);
%             else
%                 b(rowInd(ctr-1)) = b(rowInd(ctr-1)) - NA*m3ToLiter*poro(kk)*e*vale(1,jj)*xStripMiddle(rowInd(ctr-1)+jj);
%             end
%         end
%         
%         % RHS Vector entries
%         b(rowInd(ctr-1)) = b(rowInd(ctr-1)) ...
%             - rightPotentialTerm/dxC(kk)*(xStripMiddle(rowInd(ctr-1)+2*(nSpecies+1)) - xStripMiddle(rowInd(ctr-1)+(nSpecies+1))) ...
%             + leftPotentialTerm/dxC(kk)*(xStripMiddle(rowInd(ctr-1)+(nSpecies+1)) - xStripMiddle(rowInd(ctr-1)));
%         
%         b(rowInd(ctr-1)) = b(rowInd(ctr-1)) - NA*m3ToLiter*poro(kk)*e*backCharge(kk);
%     end
    
end

%------------- END OF CODE --------------