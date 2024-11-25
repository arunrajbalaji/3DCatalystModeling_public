function [dxStripConcentrationOnly, t, residual] = doTimeStep_concentrationWithHomogeneousBackup(dxC, dxF, nSpecies, xStripMiddle, xStripOld, xStripOldOld, xStripOldOldOld, isFrontNoFluxCondition, isBackNoFluxCondition,  ...
    t, dt, dtPreviousTimeStep, dtPrevPrevTimeStep, constants, rowInd, colInd, nnz, nEq, faradaicRxnInfo, isNeutralGas, faradaicRxnSourceTermsStrip, faradaicRxnSourceTermsStrip_init, transportRHSTermsStrip, transportRHSTermsStrip_init, homogeneousRxnSourceTermsStrip, homogeneousRxnSourceTermsStrip_init, stepCtr, debugMode)
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
diff = constants.diff;
acti = constants.acti;
vale = constants.vale;
e = constants.e;
kb = constants.kb;
T = constants.T;
NA = constants.NA;
phiElectrode = constants.phiElectrode;

stericACubed = constants.stericACubed;
stericOnOffVec = constants.stericOnOffVec;

xN = length(dxC);

vals = zeros(nnz,1);

t = t + dt;

% Must zero out b in advance (for reaction terms)
b = zeros(nEq,1);

% Determine LHS matrix A values and RHS vector b values
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
        
        % Unsteady term computation
        if stepCtr == 1
            timeTerm = poro(xLocIndex)/dt;
            lhsMult = 1/2;
        elseif stepCtr == 2
            timeGamma = 1/(dtPreviousTimeStep + dtPreviousTimeStep^2/dt);
            timeBeta = -timeGamma*(dt + dtPreviousTimeStep)^2/dt^2;
            timeAlpha = -timeGamma - timeBeta;
            timeTerm = poro(xLocIndex) * timeAlpha;
            lhsMult = 1;
        else
            delta_1 = dt;
            delta_2 = dt + dtPreviousTimeStep;
            delta_3 = dt + dtPreviousTimeStep + dtPrevPrevTimeStep;
            timeBeta = -delta_2*delta_3/(delta_1*(delta_1 - delta_2)*(delta_1 - delta_3));
            timeGamma = delta_3*delta_1/(delta_2*(delta_1 - delta_2)*(delta_2 - delta_3));
            timeDelta = delta_1*delta_2/(delta_3*(delta_1 - delta_3)*(delta_3 - delta_2));
            timeAlpha = -timeBeta - timeGamma - timeDelta;
            timeTerm = poro(xLocIndex) * timeAlpha;
            lhsMult = 1;
        end
        
        % species cell to left
        if xLocIndex > 2
            leftElectromigrationTerm = nextLeftElectromigrationTerm;
            leftDiffTerm = nextLeftDiffTerm;
            leftStericTerm = nextLeftStericTerm;
            vals(ctr) = lhsMult * outerDerivativeFrontFace * (leftElectromigrationTerm - leftDiffTerm*acti(xLocIndex-1,speciesIndex) - leftStericTerm);
            ctr = ctr + 1;
        end
        
        % species cell to right
        if xLocIndex < xN-1
            rightElectromigrationTerm = 1/2*((vale(1,speciesIndex)*e/kb/T)*(xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1) + (nSpecies+1)) - xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1))))/dxF(xLocIndex);
            rightDiffTerm = 1/dxF(xLocIndex)*gammaInvRight;
            rightStericTerm = stericOnOff * 1/2 * (log(1 - stericACubed * xRightStericSum) ...
                - log(1 - stericACubed * xCenterStericSum))/dxF(xLocIndex);
            vals(ctr) = lhsMult * outerDerivativeBackFace * (-rightElectromigrationTerm - rightDiffTerm*acti(xLocIndex+1,speciesIndex) + rightStericTerm);
            nextLeftElectromigrationTerm = rightElectromigrationTerm;
            nextLeftDiffTerm = rightDiffTerm;
            nextLeftStericTerm = rightStericTerm;
            
            ctr = ctr + 1;
        end
        
        % species cell center
        if xLocIndex == 2
            leftElectromigrationTerm = 1/2*((vale(1,speciesIndex)*e/kb/T)*(xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1)  + nSpecies + 1+(nSpecies+1)) - xStripMiddle(nSpecies+1)))/dxF(xLocIndex-1);
            leftDiffTerm = 1/dxF(xLocIndex-1)*gammaInvLeft;
            leftStericTerm = stericOnOff * 1/2 * (log(1 - stericACubed * xCenterStericSum) ...
                - log(1 - stericACubed * xLeftStericSum))/dxF(xLocIndex-1);
        elseif xLocIndex == xN-1
            rightElectromigrationTerm = 1/2*((vale(1,speciesIndex)*e/kb/T)*(xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1) + (nSpecies+1)) - xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1))))/dxF(xLocIndex);
            rightDiffTerm = 1/dxF(xLocIndex)*gammaInvRight;
            rightStericTerm = stericOnOff * 1/2 * (log(1 - stericACubed * xRightStericSum) ...
                - log(1 - stericACubed * xCenterStericSum))/dxF(xLocIndex);
        end
        
        if xLocIndex == 2 && isFrontNoFluxCondition && ~isNeutralGas(speciesIndex)
            vals(ctr) = timeTerm ...
                + lhsMult * outerDerivativeBackFace*(-rightElectromigrationTerm + rightDiffTerm*acti(xLocIndex,speciesIndex) + rightStericTerm);
        elseif xLocIndex == xN-1 && isBackNoFluxCondition && ~isNeutralGas(speciesIndex)
            vals(ctr) = timeTerm ...
                + lhsMult * outerDerivativeFrontFace*(leftElectromigrationTerm + leftDiffTerm*acti(xLocIndex,speciesIndex) - leftStericTerm);
        else
            vals(ctr) = timeTerm ...
                + lhsMult * outerDerivativeBackFace*(-rightElectromigrationTerm + rightDiffTerm*acti(xLocIndex,speciesIndex) + rightStericTerm)...
                + lhsMult * outerDerivativeFrontFace*(leftElectromigrationTerm + leftDiffTerm*acti(xLocIndex,speciesIndex) - leftStericTerm);
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
        
        for stericIndex=1:nSpecies
            % Cell to left
            if xLocIndex > 2
                vals(ctr) = lhsMult * lStericCrossVal;
                
                ctr = ctr + 1;
            end
            
            % Cell to right
            if xLocIndex < xN-1
                vals(ctr) = lhsMult * rStericCrossVal;
                
                ctr = ctr + 1;
            end
            
            % Cell center
            if xLocIndex == 2 && isFrontNoFluxCondition && ~isNeutralGas(speciesIndex)
                vals(ctr) = lhsMult * outerDerivativeBackFace*rightCrossSpecTermCenter;
            elseif xLocIndex == xN-1 && isBackNoFluxCondition && ~isNeutralGas(speciesIndex)
                vals(ctr) = lhsMult * outerDerivativeFrontFace*leftCrossSpecTermCenter;
            else
                vals(ctr) = lhsMult * centerStericCrossVal;
            end
            ctr = ctr + 1;
        end
    end
    
    if ~(isreal(vals) && ~any(isnan(vals)))
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
            
            vals(ctr) = lhsMult * poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
            ctr = ctr + 1;
            
            vals(ctr) = lhsMult * poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
            ctr = ctr + 1;
            
            vals(ctr) = lhsMult * poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
            ctr = ctr + 1;
            
            vals(ctr) = lhsMult * poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
            ctr = ctr + 1;
            
            if rxns(rxnIndex,3) ~= 0
                vals(ctr) = -lhsMult * poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
                ctr = ctr + 1;
                
                vals(ctr) = -lhsMult * poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
                ctr = ctr + 1;
            end
            
            if rxns(rxnIndex,4) ~= 0
                vals(ctr) = -lhsMult * poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
                ctr = ctr + 1;
                
                vals(ctr) = -lhsMult * poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
                ctr = ctr + 1;
            end
            
        elseif ((rxns(rxnIndex,1) ~= 0) && (rxns(rxnIndex,2) == 0))
            vals(ctr) = lhsMult * poro(xLocIndex)*rxns(rxnIndex,5);
            ctr = ctr + 1;
            
            if rxns(rxnIndex,3) ~= 0
                vals(ctr) = -lhsMult * poro(xLocIndex)*rxns(rxnIndex,5);
                ctr = ctr + 1;
            end
            
            if rxns(rxnIndex,4) ~= 0
                vals(ctr) = -lhsMult * poro(xLocIndex)*rxns(rxnIndex,5);
                ctr = ctr + 1;
               
            end
            
        elseif ((rxns(rxnIndex,1) == 0) && (rxns(rxnIndex,2) ~= 0))
            vals(ctr) = lhsMult * poro(xLocIndex)*rxns(rxnIndex,5);
            ctr = ctr + 1;
           
            if rxns(rxnIndex,3) ~= 0
                vals(ctr) = -lhsMult * poro(xLocIndex)*rxns(rxnIndex,5);
                ctr = ctr + 1;
            end
            
            if rxns(rxnIndex,4) ~= 0
                vals(ctr) = -lhsMult * poro(xLocIndex)*rxns(rxnIndex,5);
                ctr = ctr + 1;
            end
        end
    end
    
    % Loop over Faradaic reactions and add to right and left side as
    % necessary.
    for fRxnIndex = 1:length(faradaicRxnInfo)
        if faradaicRxnInfo(fRxnIndex).reactants ~= 0
            vals(ctr) =  0* lhsMult * 1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleReactantsPerMoleElectron*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                * 1/faradaicRxnInfo(fRxnIndex).cRef;
            ctr = ctr + 1;
        end
        
        if faradaicRxnInfo(fRxnIndex).products(1) ~= 0
            if faradaicRxnInfo(fRxnIndex).reactants ~= 0
                vals(ctr) = -0*lhsMult * 1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(1)*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                    * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                    * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                    * 1/faradaicRxnInfo(fRxnIndex).cRef;
                ctr = ctr + 1;
            end
        end
        
        if faradaicRxnInfo(fRxnIndex).products(2) ~= 0
            if faradaicRxnInfo(fRxnIndex).reactants ~= 0
                vals(ctr) = -0*lhsMult * 1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(2)*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                    * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                    * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                    * 1/faradaicRxnInfo(fRxnIndex).cRef;
                ctr = ctr + 1;
            end
        end
    end
end

% chargeField = zeros(xN-2, 1);
for sIndex = 1:nSpecies
    if stepCtr == 1
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) - poro(1)/dt*(xStripMiddle((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)) - xStripOld((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)));
    elseif stepCtr == 2
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) - poro(1)*(timeAlpha*xStripMiddle((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)) ...
            + timeBeta*xStripOld((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)) ...
            + timeGamma*xStripOldOld((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)));
    else
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) - poro(1)*(timeAlpha*xStripMiddle((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)) ...
            + timeBeta*xStripOld((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)) ...
            + timeGamma*xStripOldOld((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)) ...
            + timeDelta*xStripOldOldOld((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)));
    end

    if stepCtr == 1
%         b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) ...
%             + 1/2*reshape(faradaicRxnSourceTermsStrip(:, 1, 1, sIndex), [xN-2,1]) ...
%             + 1/2*reshape(faradaicRxnSourceTermsStrip_init(:, 1, 1, sIndex), [xN-2,1]);
        
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) ...
            + 1/2*reshape(transportRHSTermsStrip(:, 1, 1, sIndex), [xN-2,1]) ...
            + 1/2*reshape(transportRHSTermsStrip_init(:, 1, 1, sIndex), [xN-2,1]);

        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) ...
            + 1/2*reshape(homogeneousRxnSourceTermsStrip(:, 1, 1, sIndex), [xN-2,1]) ...
            + 1/2*reshape(homogeneousRxnSourceTermsStrip_init(:, 1, 1, sIndex), [xN-2,1]);
    else
%         b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(faradaicRxnSourceTermsStrip(:, 1, 1, sIndex), [xN-2,1]);
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(transportRHSTermsStrip(:, 1, 1, sIndex), [xN-2,1]);
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(homogeneousRxnSourceTermsStrip(:, 1, 1, sIndex), [xN-2,1]);
    end
    
%         chargeField = chargeField + e * vale(1, sIndex) * homogeneousRxnSourceTermsStrip(:, 1, 1, sIndex);
        
end

% chargeNorm = sqrt(sum(chargeField.^2))

%% Solve for the delta and update
% spparms('spumoni', 0)
%     warning('off', 'MATLAB:nearlySingularMatrix')
if isreal(vals) && isreal(b) && ~any(isnan(vals)) && ~any(isnan(b))
    A = sparse(rowInd, colInd, vals, nEq, nEq, length(vals));
    dxStripConcentrationOnly = A\b;
%     dxStripConcentrationOnly = bicgstab(A,b,[],1000);
%     L = ilu(A);
%     dxStripConcentrationOnly = bicg(A,b,1e-7,1000,L,L');
else
    dxStripConcentrationOnly = NaN(size(xStripMiddle));
    residual = NaN(size(xStripMiddle));
    return
end

residual = -b;

% if debugMode
%         [~, maxIndex] = max(abs(dxStripConcentrationOnly));
%         maxDelta = dxStripConcentrationOnly(maxIndex);
% else
%         maxDelta = 0;
% end

%------------- END OF CODE --------------
