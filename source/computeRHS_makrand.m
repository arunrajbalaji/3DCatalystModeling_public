function [faradaicRxnSourceTerms, implicitTransportRHSTerms, homogeneousRxnSourceTerms, explicitTransportRHSTerms] = computeRHS_makrand(dxC, dxF, dyC, dyF, dzC, dzF, dt, dtPrev, dtPrevPrev, dtPrevPrevPrev, ...
    nSpecies, fields3DStar, fields3DOld, fields3DOldOld, fields3DOldOldOld, fields3DOldOldOldOld, isFrontNoFluxCondition, isBackNoFluxCondition, isLeftNoFluxCondition, isRightNoFluxCondition, isBottomNoFluxCondition, isTopNoFluxCondition, ...
    constants, yLeftBCVec, yRightBCVec, zBottomBCVec, zTopBCVec, faradaicRxnInfo, isNeutralGas, yRightElectrolyteOverride, bulkElectrolyteL, bulkElectrolyteC, nIter, debugMode, stepCtr)
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
vale = constants.vale;
e = constants.e;
kb = constants.kb;
T = constants.T;
NA = constants.NA;
phiElectrode = constants.phiElectrode;

effectiveDiffVec = diff(1,:)' * poro(1) / tort(1);

stericACubed = constants.stericACubed;
stericOnOffVec = constants.stericOnOffVec;

xN = length(dxC);
yN = length(dyC);
zN = length(dzC);

nVar = nSpecies+1;

fields3DStar = cat(2, yLeftBCVec, fields3DStar, yRightBCVec);
fields3DStar = cat(3, cat(2, zBottomBCVec(:, 1, 1), zBottomBCVec, zBottomBCVec(:, end, 1)), ...
    fields3DStar, ...
    cat(2, zTopBCVec(:, 1, end), zTopBCVec, zTopBCVec(:, end, end)));

fields3DOld = cat(2, yLeftBCVec, fields3DOld, yRightBCVec);
fields3DOld = cat(3, cat(2, zBottomBCVec(:, 1, 1), zBottomBCVec, zBottomBCVec(:, end, 1)), ...
    fields3DOld, ...
    cat(2, zTopBCVec(:, 1, end), zTopBCVec, zTopBCVec(:, end, end)));

% fields3DOldOld = cat(2, yLeftBCVec, fields3DOldOld, yRightBCVec);
% fields3DOldOld = cat(3, cat(2, zBottomBCVec(:, 1, 1), zBottomBCVec, zBottomBCVec(:, end, 1)), ...
%     fields3DOldOld, ...
%     cat(2, zTopBCVec(:, 1, end), zTopBCVec, zTopBCVec(:, end, end)));
% 
% fields3DOldOldOld = cat(2, yLeftBCVec, fields3DOldOldOld, yRightBCVec);
% fields3DOldOldOld = cat(3, cat(2, zBottomBCVec(:, 1, 1), zBottomBCVec, zBottomBCVec(:, end, 1)), ...
%     fields3DOldOldOld, ...
%     cat(2, zTopBCVec(:, 1, end), zTopBCVec, zTopBCVec(:, end, end)));
% 
% fields3DOldOldOldOld = cat(2, yLeftBCVec, fields3DOldOldOldOld, yRightBCVec);
% fields3DOldOldOldOld = cat(3, cat(2, zBottomBCVec(:, 1, 1), zBottomBCVec, zBottomBCVec(:, end, 1)), ...
%     fields3DOldOldOldOld, ...
%     cat(2, zTopBCVec(:, 1, end), zTopBCVec, zTopBCVec(:, end, end)));

% if stepCtr == 1
%     timeAlpha = 1;
%     timeBeta = 0;
%     timeGamma = 0;
%     timeDelta = 0;
% elseif stepCtr == 2
%     timeBeta = -dt/dtPrev;
%     timeAlpha = (dt + dtPrev)/dtPrev;
%     timeGamma = 0;
%     timeDelta = 0;
% elseif stepCtr == 3
%     timeBeta = -((dt + dtPrev + dtPrevPrev)^2 + (dt + dtPrev + dtPrevPrev)/(dtPrev + dtPrevPrev) * (dt^2 - (dt + dtPrev + dtPrevPrev)^2))/((dt + dtPrev)^2 - (dt + dtPrev + dtPrevPrev)^2 - dtPrevPrev/(dtPrev + dtPrevPrev)*(dt^2 - (dt + dtPrev + dtPrevPrev)^2));
%     timeAlpha = (dt + dtPrev + dtPrevPrev)/(dtPrev + dtPrevPrev) - dtPrevPrev/(dtPrev + dtPrevPrev) * timeBeta;
%     timeGamma = 1 - timeAlpha - timeBeta;
%     timeDelta = 0;
% else
%    timeLHSVec = [dt (dt+dtPrev) (dt+dtPrev+dtPrevPrev) (dt+dtPrev+dtPrevPrev+dtPrevPrevPrev)];
%    timeLHS = [1 1 1 1; timeLHSVec; timeLHSVec.^2; timeLHSVec.^3];
%    timeRHS = [1; 0; 0; 0];
%    timeSolve = timeLHS \ timeRHS;
%    timeAlpha = timeSolve(1);
%    timeBeta = timeSolve(2);
%    timeGamma = timeSolve(3);
%    timeDelta = timeSolve(4);
% end

faradaicRxnSourceTerms = zeros(xN-2, yN, zN, nSpecies);
implicitTransportRHSTerms = zeros(xN-2, yN, zN, nSpecies);
homogeneousRxnSourceTerms = zeros(xN-2, yN, zN, nSpecies);
explicitTransportRHSTerms = zeros(xN-2, yN, zN, nSpecies);

stericSum = zeros(xN, yN+2, zN+2);

% dPhidX = (fields3DStar(2*nVar:nVar:end, 2:end-1, 2:end-1) - fields3DStar(nVar:nVar:end-nVar, 2:end-1, 2:end-1))./dxF;

% k1Frac = rhsFunctionTransverse(fields3DOld, dxC, dyC, dyF, dzC, dzF, constants, isLeftNoFluxCondition, isRightNoFluxCondition, isBottomNoFluxCondition, isTopNoFluxCondition, yRightElectrolyteOverride, bulkElectrolyteL, bulkElectrolyteC, isNeutralGas, nSpecies);
% 
% k2Frac = rhsFunctionTransverse(fields3DOld + dt * 1/2*k1Frac, dxC, dyC, dyF, dzC, dzF, constants, isLeftNoFluxCondition, isRightNoFluxCondition, isBottomNoFluxCondition, isTopNoFluxCondition, yRightElectrolyteOverride, bulkElectrolyteL, bulkElectrolyteC, isNeutralGas, nSpecies);
% 
% k3Frac = rhsFunctionTransverse(fields3DOld + dt * 1/2*k2Frac, dxC, dyC, dyF, dzC, dzF, constants, isLeftNoFluxCondition, isRightNoFluxCondition, isBottomNoFluxCondition, isTopNoFluxCondition, yRightElectrolyteOverride, bulkElectrolyteL, bulkElectrolyteC, isNeutralGas, nSpecies);
% 
% k4Frac = rhsFunctionTransverse(fields3DOld + dt * k3Frac, dxC, dyC, dyF, dzC, dzF, constants, isLeftNoFluxCondition, isRightNoFluxCondition, isBottomNoFluxCondition, isTopNoFluxCondition, yRightElectrolyteOverride, bulkElectrolyteL, bulkElectrolyteC, isNeutralGas, nSpecies);

explicitRHSAllSpeciesOld = rhsFunctionTransverse(fields3DOld, dxC, dyC, dyF, dzC, dzF, constants, isLeftNoFluxCondition, isRightNoFluxCondition, isBottomNoFluxCondition, isTopNoFluxCondition, yRightElectrolyteOverride, bulkElectrolyteL, bulkElectrolyteC, isNeutralGas, nSpecies);
explicitRHSAllSpeciesOldOld = rhsFunctionTransverse(fields3DOldOld, dxC, dyC, dyF, dzC, dzF, constants, isLeftNoFluxCondition, isRightNoFluxCondition, isBottomNoFluxCondition, isTopNoFluxCondition, yRightElectrolyteOverride, bulkElectrolyteL, bulkElectrolyteC, isNeutralGas, nSpecies);

for xIndex = 1:xN
    stericSum(xIndex, :, :) = sum(fields3DStar((xIndex-1)*nVar+1:(xIndex-1)*nVar+nSpecies, :, :) .* stericOnOffVec, 1);
end

for sIndex = 1:nSpecies
    % Compute right hand side - X Fluxes
    xDiff = -effectiveDiffVec(sIndex) * (fields3DStar(nVar+sIndex:nVar:end, 2:end-1, 2:end-1) - fields3DStar(sIndex:nVar:end-nVar, 2:end-1, 2:end-1)) ./ dxF;
    xElectro = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * (fields3DStar(nVar+sIndex:nVar:end, 2:end-1, 2:end-1) + fields3DStar(sIndex:nVar:end-nVar, 2:end-1, 2:end-1))/2 ...
        .* (fields3DStar(2*nVar:nVar:end, 2:end-1, 2:end-1) - fields3DStar(nVar:nVar:end-nVar, 2:end-1, 2:end-1)) ./ dxF;
    xSteric = -0*effectiveDiffVec(sIndex) * (fields3DStar(nVar+sIndex:nVar:end, 2:end-1, 2:end-1) + fields3DStar(sIndex:nVar:end-nVar, 2:end-1, 2:end-1))/2 ...
        .* (log(1 - stericACubed * stericSum(2:end, 2:end-1, 2:end-1)) - log(1 - stericACubed * stericSum(1:end-1, 2:end-1, 2:end-1))) ./ dxF;

    xDiffOld = -effectiveDiffVec(sIndex) * (fields3DOld(nVar+sIndex:nVar:end, 2:end-1, 2:end-1) - fields3DOld(sIndex:nVar:end-nVar, 2:end-1, 2:end-1)) ./ dxF;
    xElectroOld = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * (fields3DOld(nVar+sIndex:nVar:end, 2:end-1, 2:end-1) + fields3DOld(sIndex:nVar:end-nVar, 2:end-1, 2:end-1))/2 ...
        .* (fields3DOld(2*nVar:nVar:end, 2:end-1, 2:end-1) - fields3DOld(nVar:nVar:end-nVar, 2:end-1, 2:end-1)) ./ dxF;
    xStericOld = -0*effectiveDiffVec(sIndex) * (fields3DOld(nVar+sIndex:nVar:end, 2:end-1, 2:end-1) + fields3DOld(sIndex:nVar:end-nVar, 2:end-1, 2:end-1))/2 ...
        .* (log(1 - stericACubed * stericSum(2:end, 2:end-1, 2:end-1)) - log(1 - stericACubed * stericSum(1:end-1, 2:end-1, 2:end-1))) ./ dxF;

%     emOperator_1 = diag(1/2*(dPhidX(1:end-1, 5, 5) + dPhidX(2:end, 5, 5))) ...
%         * (diag(1./(dxF(1:end-2)+dxF(2:end-1)), 1) + diag(-1./(dxF(2:end-1)+dxF(3:end)), -1));
% 
%     emOperator_2 = diag(((fields3DStar(3*nVar:nVar:end, 5, 5) - fields3DStar(2*nVar:nVar:end-nVar, 5, 5))./dxF(2:end) ...
%         - (fields3DStar(2*nVar:nVar:end-nVar, 5, 5) - fields3DStar(nVar:nVar:end-2*nVar, 5, 5))./dxF(1:end-1)) ./ dxC(2:end-1));
% 
%     totalEmOperator = diff(1, sIndex) * vale(1,sIndex) * e / (kb * T) * (emOperator_1 + emOperator_2);

%     xElectroDelta = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * (fields3DStar(nVar+sIndex:nVar:end, 2:end-1, 2:end-1) + fields3DStar(sIndex:nVar:end-nVar, 2:end-1, 2:end-1))/2 ...
%         .* (deltaFields(2*nVar:nVar:end, 2:end-1, 2:end-1) - deltaFields(nVar:nVar:end-nVar, 2:end-1, 2:end-1)) ./ dxF;
    
    xTotalFluxStar = xDiff + xElectro + xSteric;% + xElectroDelta;
    xTotalFluxOld = xDiffOld + xElectroOld + xStericOld;% + xElectroDelta;

    if isFrontNoFluxCondition && ~isNeutralGas(sIndex)
        xTotalFluxStar(1, :, :) = 0;
        xTotalFluxOld(1, :, :) = 0;
    end

    if isBackNoFluxCondition && ~isNeutralGas(sIndex)
        xTotalFluxStar(end,: , :) = 0;
        xTotalFluxOld(end,: , :) = 0;
    end

    % Override for neutral gasses (correct Dirichlet BC)
    if isNeutralGas(sIndex)
        xTotalFluxStar(1, :, :) = 2 * effectiveDiffVec(sIndex)/dxF(1) * (fields3DStar(sIndex, 2:end-1, 2:end-1) ...
            - fields3DStar(nVar+sIndex, 2:end-1, 2:end-1));
        xTotalFluxStar(end, :, :) = 2 * effectiveDiffVec(sIndex)/dxF(end) *(-fields3DStar(end-nVar+sIndex, 2:end-1, 2:end-1) ...
            + fields3DStar(end-2*nVar+sIndex, 2:end-1, 2:end-1));


        xTotalFluxOld(1, :, :) = 2 * effectiveDiffVec(sIndex)/dxF(1) * (fields3DOld(sIndex, 2:end-1, 2:end-1) ...
            - fields3DOld(nVar+sIndex, 2:end-1, 2:end-1));
        xTotalFluxOld(end, :, :) = 2 * effectiveDiffVec(sIndex)/dxF(end) *(-fields3DOld(end-nVar+sIndex, 2:end-1, 2:end-1) ...
            + fields3DOld(end-2*nVar+sIndex, 2:end-1, 2:end-1));

%         xTotalFlux(1,: , :) = 0;
%         xTotalFlux(end,: , :) = 0;
    end

%     if yRightElectrolyteOverride
%         if (vale(1, sIndex) > 0) && (sIndex == 3)   % HYDROGEN
%             K_eq = constants.rxns{1}(15,5)/constants.rxns{1}(16,5);
%             % Find OH- concentration at ghost cell center
%             ghostCellOH = -fields3DOld(7:nVar:end, end-2, :) + 2 * fields3DOld(7:nVar:end, end-1, :);
%             % Equilibrate water recombinationm, Set Dirichlet value
%             fields3DOld(7:nVar:end, end, :) = K_eq ./ ghostCellOH;
%             
%             % Find OH- concentration at ghost cell center
%             ghostCellOH = -fields3DOldOld(7:nVar:end, end-2, :) + 2 * fields3DOldOld(7:nVar:end, end-1, :);
%             % Equilibrate water recombinationm, Set Dirichlet value
%             fields3DOldOld(7:nVar:end, end, :) = K_eq ./ ghostCellOH;
%             
%             % Find OH- concentration at ghost cell center
%             ghostCellOH = -fields3DOldOldOld(7:nVar:end, end-2, :) + 2 * fields3DOldOldOld(7:nVar:end, end-1, :);
%             % Equilibrate water recombinationm, Set Dirichlet value
%             fields3DOldOldOld(7:nVar:end, end, :) = K_eq ./ ghostCellOH;
%             
%             % Find OH- concentration at ghost cell center
%             ghostCellOH = -fields3DOldOldOldOld(7:nVar:end, end-2, :) + 2 * fields3DOldOldOldOld(7:nVar:end, end-1, :);
%             % Equilibrate water recombinationm, Set Dirichlet value
%             fields3DOldOldOldOld(7:nVar:end, end, :) = K_eq ./ ghostCellOH;
%         elseif (vale(1, sIndex) > 0) && (sIndex == 6)   % POTASSIUM
%             % DO NOTHING (just use Dirichlet Eq. solution value)
%         end
%     end
% 
%     % Compute right hand side - Y Fluxes
%     yDiffOld = -effectiveDiffVec(sIndex) * (fields3DOld(nVar+sIndex:nVar:end-nVar, 2:end, 2:end-1) - fields3DOld(nVar+sIndex:nVar:end-nVar, 1:end-1, 2:end-1)) ./ reshape(dyF, [1, yN+1]);
%     yElectroOld = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * (fields3DOld(nVar+sIndex:nVar:end-nVar, 2:end, 2:end-1) + fields3DOld(nVar+sIndex:nVar:end-nVar, 1:end-1, 2:end-1))/2 ...
%         .* (fields3DOld(2*nVar:nVar:end-nVar, 2:end, 2:end-1) - (fields3DOld(2*nVar:nVar:end-nVar, 1:end-1, 2:end-1))) ./ reshape(dyF, [1, yN+1]);
%     % Correct Dirichlet treatment for yRight boundary
%     yDiffOld(:, end, :) = -effectiveDiffVec(sIndex) * 2 * (fields3DOld(nVar+sIndex:nVar:end-nVar, end, 2:end-1) - fields3DOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1)) ./ dyF(end);
%     yElectroOld(:, end, :) = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * fields3DOld(nVar+sIndex:nVar:end-nVar, end, 2:end-1) ...
%         * 2 .* (fields3DOld(2*nVar:nVar:end-nVar, end, 2:end-1) - fields3DOld(2*nVar:nVar:end-nVar, end-1, 2:end-1)) ./ dyF(end);
%     yStericOld = -0*effectiveDiffVec(sIndex) * (fields3DOld(nVar+sIndex:nVar:end-nVar, 2:end, 2:end-1) + fields3DOld(nVar+sIndex:nVar:end-nVar, 1:end-1, 2:end-1))/2 ...
%         .* (log(1 - stericACubed * stericSum(2:end-1, 2:end, 2:end-1)) - log(1 - stericACubed * stericSum(2:end-1, 1:end-1, 2:end-1))) ./ reshape(dyF, [1, yN+1]);
%     
%     % Compute right hand side - Y Fluxes
%     yDiffOldOld = -effectiveDiffVec(sIndex) * (fields3DOldOld(nVar+sIndex:nVar:end-nVar, 2:end, 2:end-1) - fields3DOldOld(nVar+sIndex:nVar:end-nVar, 1:end-1, 2:end-1)) ./ reshape(dyF, [1, yN+1]);
%     yElectroOldOld = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * (fields3DOldOld(nVar+sIndex:nVar:end-nVar, 2:end, 2:end-1) + fields3DOldOld(nVar+sIndex:nVar:end-nVar, 1:end-1, 2:end-1))/2 ...
%         .* (fields3DOldOld(2*nVar:nVar:end-nVar, 2:end, 2:end-1) - (fields3DOldOld(2*nVar:nVar:end-nVar, 1:end-1, 2:end-1))) ./ reshape(dyF, [1, yN+1]);
%     % Correct Dirichlet treatment for yRight boundary
%     yDiffOldOld(:, end, :) = -effectiveDiffVec(sIndex) * 2 * (fields3DOldOld(nVar+sIndex:nVar:end-nVar, end, 2:end-1) - fields3DOldOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1)) ./ dyF(end);
%     yElectroOldOld(:, end, :) = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * fields3DOldOld(nVar+sIndex:nVar:end-nVar, end, 2:end-1) ...
%         * 2 .* (fields3DOldOld(2*nVar:nVar:end-nVar, end, 2:end-1) - fields3DOldOld(2*nVar:nVar:end-nVar, end-1, 2:end-1)) ./ dyF(end);
%     yStericOldOld = -0*effectiveDiffVec(sIndex) * (fields3DOldOld(nVar+sIndex:nVar:end-nVar, 2:end, 2:end-1) + fields3DOldOld(nVar+sIndex:nVar:end-nVar, 1:end-1, 2:end-1))/2 ...
%         .* (log(1 - stericACubed * stericSum(2:end-1, 2:end, 2:end-1)) - log(1 - stericACubed * stericSum(2:end-1, 1:end-1, 2:end-1))) ./ reshape(dyF, [1, yN+1]);
%     
%     % Compute right hand side - Y Fluxes
%     yDiffOldOldOld = -effectiveDiffVec(sIndex) * (fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end, 2:end-1) - fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, 1:end-1, 2:end-1)) ./ reshape(dyF, [1, yN+1]);
%     yElectroOldOldOld = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * (fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end, 2:end-1) + fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, 1:end-1, 2:end-1))/2 ...
%         .* (fields3DOldOldOld(2*nVar:nVar:end-nVar, 2:end, 2:end-1) - (fields3DOldOldOld(2*nVar:nVar:end-nVar, 1:end-1, 2:end-1))) ./ reshape(dyF, [1, yN+1]);
%     % Correct Dirichlet treatment for yRight boundary
%     yDiffOldOldOld(:, end, :) = -effectiveDiffVec(sIndex) * 2 * (fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, end, 2:end-1) - fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1)) ./ dyF(end);
%     yElectroOldOldOld(:, end, :) = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, end, 2:end-1) ...
%         * 2 .* (fields3DOldOldOld(2*nVar:nVar:end-nVar, end, 2:end-1) - fields3DOldOldOld(2*nVar:nVar:end-nVar, end-1, 2:end-1)) ./ dyF(end);
%     yStericOldOldOld = -0*effectiveDiffVec(sIndex) * (fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end, 2:end-1) + fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, 1:end-1, 2:end-1))/2 ...
%         .* (log(1 - stericACubed * stericSum(2:end-1, 2:end, 2:end-1)) - log(1 - stericACubed * stericSum(2:end-1, 1:end-1, 2:end-1))) ./ reshape(dyF, [1, yN+1]);
%     
%     % Compute right hand side - Y Fluxes
%     yDiffOldOldOldOld = -effectiveDiffVec(sIndex) * (fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end, 2:end-1) - fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, 1:end-1, 2:end-1)) ./ reshape(dyF, [1, yN+1]);
%     yElectroOldOldOldOld = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * (fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end, 2:end-1) + fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, 1:end-1, 2:end-1))/2 ...
%         .* (fields3DOldOldOldOld(2*nVar:nVar:end-nVar, 2:end, 2:end-1) - (fields3DOldOldOldOld(2*nVar:nVar:end-nVar, 1:end-1, 2:end-1))) ./ reshape(dyF, [1, yN+1]);
%     % Correct Dirichlet treatment for yRight boundary
%     yDiffOldOldOldOld(:, end, :) = -effectiveDiffVec(sIndex) * 2 * (fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, end, 2:end-1) - fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1)) ./ dyF(end);
%     yElectroOldOldOldOld(:, end, :) = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, end, 2:end-1) ...
%         * 2 .* (fields3DOldOldOldOld(2*nVar:nVar:end-nVar, end, 2:end-1) - fields3DOldOldOldOld(2*nVar:nVar:end-nVar, end-1, 2:end-1)) ./ dyF(end);
%     yStericOldOldOldOld = -0*effectiveDiffVec(sIndex) * (fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end, 2:end-1) + fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, 1:end-1, 2:end-1))/2 ...
%         .* (log(1 - stericACubed * stericSum(2:end-1, 2:end, 2:end-1)) - log(1 - stericACubed * stericSum(2:end-1, 1:end-1, 2:end-1))) ./ reshape(dyF, [1, yN+1]);
% 
%     if yRightElectrolyteOverride
%         if vale(1, sIndex) < 0
%             delta_1_2 = dyC(end)/2;
%             delta_3_2 = dyC(end)/2 + dyF(end-1);
%             delta_5_2 = dyC(end)/2 + dyF(end-1) + dyF(end-2);
%             beta = (-1/(delta_5_2 - delta_1_2)*(delta_1_2^2 - delta_5_2^2)) / (-(delta_5_2 - delta_3_2)/(delta_5_2 - delta_1_2) * (delta_1_2^2 - delta_5_2^2) + delta_3_2^2 - delta_5_2^2);
%             alpha = 1/(delta_5_2 - delta_1_2) - beta * (delta_5_2 - delta_3_2)/(delta_5_2 - delta_1_2);
%             gamma = -alpha - beta;
%             
%             yDiffOld(:, end, :) = - effectiveDiffVec(sIndex) * (alpha * fields3DOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1) + beta * fields3DOld(nVar+sIndex:nVar:end-nVar, end-2, 2:end-1) + gamma * fields3DOld(nVar+sIndex:nVar:end-nVar, end-3, 2:end-1));
%             yElectroOld(:, end, :) = -vale(1,sIndex)*e/(kb*T)*effectiveDiffVec(sIndex) * (-1/2 * fields3DOld(nVar+sIndex:nVar:end-nVar, end-2, 2:end-1) + 3/2 * fields3DOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1)) ...
%                 * 2 .* (fields3DOld(2*nVar:nVar:end-nVar, end, 2:end-1) - (fields3DOld(2*nVar:nVar:end-nVar, end-1, 2:end-1) ))/dyF(end);
%             effectiveBCConc = -1 * fields3DOld(:, end-2, :) + 2 * fields3DOld(:, end-1, :);
%             fields3DOld(sIndex:nVar:end, end, :) = effectiveBCConc(sIndex:nVar:end, 1, :);
%             
%             yDiffOldOld(:, end, :) = - effectiveDiffVec(sIndex) * (alpha * fields3DOldOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1) + beta * fields3DOldOld(nVar+sIndex:nVar:end-nVar, end-2, 2:end-1) + gamma * fields3DOldOld(nVar+sIndex:nVar:end-nVar, end-3, 2:end-1));
%             yElectroOldOld(:, end, :) = -vale(1,sIndex)*e/(kb*T)*effectiveDiffVec(sIndex) * (-1/2 * fields3DOldOld(nVar+sIndex:nVar:end-nVar, end-2, 2:end-1) + 3/2 * fields3DOldOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1)) ...
%                 * 2 .* (fields3DOldOld(2*nVar:nVar:end-nVar, end, 2:end-1) - (fields3DOldOld(2*nVar:nVar:end-nVar, end-1, 2:end-1) ))/dyF(end);
%             effectiveBCConc = -1 * fields3DOldOld(:, end-2, :) + 2 * fields3DOldOld(:, end-1, :);
%             fields3DOldOld(sIndex:nVar:end, end, :) = effectiveBCConc(sIndex:nVar:end, 1, :);
%             
%             yDiffOldOldOld(:, end, :) = - effectiveDiffVec(sIndex) * (alpha * fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1) + beta * fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, end-2, 2:end-1) + gamma * fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, end-3, 2:end-1));
%             yElectroOldOldOld(:, end, :) = -vale(1,sIndex)*e/(kb*T)*effectiveDiffVec(sIndex) * (-1/2 * fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, end-2, 2:end-1) + 3/2 * fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1)) ...
%                 * 2 .* (fields3DOldOldOld(2*nVar:nVar:end-nVar, end, 2:end-1) - (fields3DOldOldOld(2*nVar:nVar:end-nVar, end-1, 2:end-1) ))/dyF(end);
%             effectiveBCConc = -1 * fields3DOldOldOld(:, end-2, :) + 2 * fields3DOldOldOld(:, end-1, :);
%             fields3DOldOldOld(sIndex:nVar:end, end, :) = effectiveBCConc(sIndex:nVar:end, 1, :);
%             
%             yDiffOldOldOldOld(:, end, :) = - effectiveDiffVec(sIndex) * (alpha * fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1) + beta * fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, end-2, 2:end-1) + gamma * fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, end-3, 2:end-1));
%             yElectroOldOldOldOld(:, end, :) = -vale(1,sIndex)*e/(kb*T)*effectiveDiffVec(sIndex) * (-1/2 * fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, end-2, 2:end-1) + 3/2 * fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1)) ...
%                 * 2 .* (fields3DOldOldOldOld(2*nVar:nVar:end-nVar, end, 2:end-1) - (fields3DOldOldOldOld(2*nVar:nVar:end-nVar, end-1, 2:end-1) ))/dyF(end);
%             effectiveBCConc = -1 * fields3DOldOldOldOld(:, end-2, :) + 2 * fields3DOldOldOldOld(:, end-1, :);
%             fields3DOldOldOldOld(sIndex:nVar:end, end, :) = effectiveBCConc(sIndex:nVar:end, 1, :);
%         elseif vale(1, sIndex) == 0
%             yDiffOld(:, end, :) = -poro(1)/tort(1)*(diff(1,1)+diff(1,4))/2 / bulkElectrolyteL * (bulkElectrolyteC(sIndex) - (-1/2 * fields3DOld(nVar+sIndex:nVar:end-nVar, end-2, 2:end-1) + 3/2 * fields3DOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1)));
%             yElectroOld(:, end, :) = 0;
%             
%             yDiffOldOld(:, end, :) = -poro(1)/tort(1)*(diff(1,1)+diff(1,4))/2 / bulkElectrolyteL * (bulkElectrolyteC(sIndex) - (-1/2 * fields3DOldOld(nVar+sIndex:nVar:end-nVar, end-2, 2:end-1) + 3/2 * fields3DOldOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1)));
%             yElectroOldOld(:, end, :) = 0;
%             
%             yDiffOldOldOld(:, end, :) = -poro(1)/tort(1)*(diff(1,1)+diff(1,4))/2 / bulkElectrolyteL * (bulkElectrolyteC(sIndex) - (-1/2 * fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, end-2, 2:end-1) + 3/2 * fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1)));
%             yElectroOldOldOld(:, end, :) = 0;
%             
%             yDiffOldOldOldOld(:, end, :) = -poro(1)/tort(1)*(diff(1,1)+diff(1,4))/2 / bulkElectrolyteL * (bulkElectrolyteC(sIndex) - (-1/2 * fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, end-2, 2:end-1) + 3/2 * fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1)));
%             yElectroOldOldOldOld(:, end, :) = 0;
% 
%         elseif vale(1, sIndex) > 0
%             yElectroOld(:, end, :) = -vale(1,sIndex)*e/(kb*T)*effectiveDiffVec(sIndex) * (fields3DOld(nVar+sIndex:nVar:end-nVar, end, 2:end-1) + fields3DOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1))/2 ...
%                 * 2 .* (fields3DOld(2*nVar:nVar:end-nVar, end, 2:end-1) - (fields3DOld(2*nVar:nVar:end-nVar, end-1, 2:end-1)))/dyF(end);
%             
%             yElectroOldOld(:, end, :) = -vale(1,sIndex)*e/(kb*T)*effectiveDiffVec(sIndex) * (fields3DOldOld(nVar+sIndex:nVar:end-nVar, end, 2:end-1) + fields3DOldOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1))/2 ...
%                 * 2 .* (fields3DOldOld(2*nVar:nVar:end-nVar, end, 2:end-1) - (fields3DOldOld(2*nVar:nVar:end-nVar, end-1, 2:end-1)))/dyF(end);
%             
%             yElectroOldOldOld(:, end, :) = -vale(1,sIndex)*e/(kb*T)*effectiveDiffVec(sIndex) * (fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, end, 2:end-1) + fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1))/2 ...
%                 * 2 .* (fields3DOldOldOld(2*nVar:nVar:end-nVar, end, 2:end-1) - (fields3DOldOldOld(2*nVar:nVar:end-nVar, end-1, 2:end-1)))/dyF(end);
%             
%             yElectroOldOldOldOld(:, end, :) = -vale(1,sIndex)*e/(kb*T)*effectiveDiffVec(sIndex) * (fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, end, 2:end-1) + fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1))/2 ...
%                 * 2 .* (fields3DOldOldOldOld(2*nVar:nVar:end-nVar, end, 2:end-1) - (fields3DOldOldOldOld(2*nVar:nVar:end-nVar, end-1, 2:end-1)))/dyF(end);
%         end
%     end
% 
%     yTotalFluxOld = yDiffOld + yElectroOld + yStericOld;
%     yTotalFluxOldOld = yDiffOldOld + yElectroOldOld + yStericOldOld;
%     yTotalFluxOldOldOld = yDiffOldOldOld + yElectroOldOldOld + yStericOldOldOld;
%     yTotalFluxOldOldOldOld = yDiffOldOldOldOld + yElectroOldOldOldOld + yStericOldOldOldOld;
% 
%     if isLeftNoFluxCondition && ~isNeutralGas(sIndex)
%         yTotalFluxOld(:, 1, :) = 0;
%         yTotalFluxOldOld(:, 1, :) = 0;
%         yTotalFluxOldOldOld(:, 1, :) = 0;
%         yTotalFluxOldOldOldOld(:, 1, :) = 0;
%     end
% 
%     if isRightNoFluxCondition && ~isNeutralGas(sIndex) && ~yRightElectrolyteOverride
%         yTotalFluxOld(:, end, :) = 0;
%         yTotalFluxOldOld(:, end, :) = 0;
%         yTotalFluxOldOldOld(:, end, :) = 0;
%         yTotalFluxOldOldOldOld(:, end, :) = 0;
%     end
% 
%     % Override for neutral gasses (correct Dirichlet BC)
%     if isNeutralGas(sIndex)
% %         yTotalFluxOld(:, 1, :) = 2 * effectiveDiffVec(sIndex)/dyF(1) * fields3DOld(nVar+sIndex:nVar:end-nVar, 1, 2:end-1) ...
% %             - 2 * effectiveDiffVec(sIndex)/dyF(1) * fields3DOld(nVar+sIndex:nVar:end-nVar, 2, 2:end-1);
%         yTotalFluxOld(:, 1, :) = 0;
%         yTotalFluxOldOld(:, 1, :) = 0;
%         yTotalFluxOldOldOld(:, 1, :) = 0;
%         yTotalFluxOldOldOldOld(:, 1, :) = 0;
%     end
% 
%     % Compute right hand side - Z Fluxes
%     zDiffOld = -effectiveDiffVec(sIndex) * (fields3DOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end) - fields3DOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 1:end-1)) ./ reshape(dzF, [1, 1, zN+1]);
%     zElectroOld = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * (fields3DOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end) + fields3DOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 1:end-1))/2 ...
%         .* (fields3DOld(2*nVar:nVar:end-nVar, 2:end-1, 2:end) - (fields3DOld(2*nVar:nVar:end-nVar, 2:end-1, 1:end-1))) ./ reshape(dzF, [1, 1, zN+1]);
%     zStericOld = -0*effectiveDiffVec(sIndex) * (fields3DOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end) + fields3DOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 1:end-1))/2 ...
%         .* (log(1 - stericACubed * stericSum(2:end-1, 2:end-1, 2:end)) - log(1 - stericACubed * stericSum(2:end-1, 2:end-1, 1:end-1))) ./ reshape(dzF, [1, 1, zN+1]);
% 
% 
%     zTotalFluxOld = zDiffOld + zElectroOld + zStericOld;
% 
%     if isBottomNoFluxCondition
%         zTotalFluxOld(:, :, 1) = 0;
%     end
% 
%     if isTopNoFluxCondition
%         zTotalFluxOld(:, :, end) = 0;
%     end
%     
%     % Compute right hand side - Z Fluxes
%     zDiffOldOld = -effectiveDiffVec(sIndex) * (fields3DOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end) - fields3DOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 1:end-1)) ./ reshape(dzF, [1, 1, zN+1]);
%     zElectroOldOld = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * (fields3DOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end) + fields3DOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 1:end-1))/2 ...
%         .* (fields3DOldOld(2*nVar:nVar:end-nVar, 2:end-1, 2:end) - (fields3DOldOld(2*nVar:nVar:end-nVar, 2:end-1, 1:end-1))) ./ reshape(dzF, [1, 1, zN+1]);
%     zStericOldOld = -0*effectiveDiffVec(sIndex) * (fields3DOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end) + fields3DOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 1:end-1))/2 ...
%         .* (log(1 - stericACubed * stericSum(2:end-1, 2:end-1, 2:end)) - log(1 - stericACubed * stericSum(2:end-1, 2:end-1, 1:end-1))) ./ reshape(dzF, [1, 1, zN+1]);
% 
% 
%     zTotalFluxOldOld = zDiffOldOld + zElectroOldOld + zStericOldOld;
% 
%     if isBottomNoFluxCondition
%         zTotalFluxOldOld(:, :, 1) = 0;
%     end
% 
%     if isTopNoFluxCondition
%         zTotalFluxOldOld(:, :, end) = 0;
%     end
%     
%     % Compute right hand side - Z Fluxes
%     zDiffOldOldOld = -effectiveDiffVec(sIndex) * (fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end) - fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 1:end-1)) ./ reshape(dzF, [1, 1, zN+1]);
%     zElectroOldOldOld = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * (fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end) + fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 1:end-1))/2 ...
%         .* (fields3DOldOldOld(2*nVar:nVar:end-nVar, 2:end-1, 2:end) - (fields3DOldOldOld(2*nVar:nVar:end-nVar, 2:end-1, 1:end-1))) ./ reshape(dzF, [1, 1, zN+1]);
%     zStericOldOldOld = -0*effectiveDiffVec(sIndex) * (fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end) + fields3DOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 1:end-1))/2 ...
%         .* (log(1 - stericACubed * stericSum(2:end-1, 2:end-1, 2:end)) - log(1 - stericACubed * stericSum(2:end-1, 2:end-1, 1:end-1))) ./ reshape(dzF, [1, 1, zN+1]);
% 
% 
%     zTotalFluxOldOldOld = zDiffOldOldOld + zElectroOldOldOld + zStericOldOldOld;
% 
%     if isBottomNoFluxCondition
%         zTotalFluxOldOldOld(:, :, 1) = 0;
%     end
% 
%     if isTopNoFluxCondition
%         zTotalFluxOldOldOld(:, :, end) = 0;
%     end
%     
%     % Compute right hand side - Z Fluxes
%     zDiffOldOldOldOld = -effectiveDiffVec(sIndex) * (fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end) - fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 1:end-1)) ./ reshape(dzF, [1, 1, zN+1]);
%     zElectroOldOldOldOld = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * (fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end) + fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 1:end-1))/2 ...
%         .* (fields3DOldOldOldOld(2*nVar:nVar:end-nVar, 2:end-1, 2:end) - (fields3DOldOldOldOld(2*nVar:nVar:end-nVar, 2:end-1, 1:end-1))) ./ reshape(dzF, [1, 1, zN+1]);
%     zStericOldOldOldOld = -0*effectiveDiffVec(sIndex) * (fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end) + fields3DOldOldOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 1:end-1))/2 ...
%         .* (log(1 - stericACubed * stericSum(2:end-1, 2:end-1, 2:end)) - log(1 - stericACubed * stericSum(2:end-1, 2:end-1, 1:end-1))) ./ reshape(dzF, [1, 1, zN+1]);
% 
% 
%     zTotalFluxOldOldOldOld = zDiffOldOldOldOld + zElectroOldOldOldOld + zStericOldOldOldOld;
% 
%     if isBottomNoFluxCondition
%         zTotalFluxOldOldOld(:, :, 1) = 0;
%     end
% 
%     if isTopNoFluxCondition
%         zTotalFluxOldOldOld(:, :, end) = 0;
%     end

    if stepCtr == 1
        explicitTransportRHSTerms(:, :, :, sIndex) = explicitRHSAllSpeciesOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end-1);
    else
        explicitTransportRHSTerms(:, :, :, sIndex) = ...
            (dt + dtPrev)/dtPrev * explicitRHSAllSpeciesOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end-1) ...
            - dt/dtPrev * explicitRHSAllSpeciesOldOld(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end-1);
    end
    
    if stepCtr == 1
        implicitTransportRHSTerms(:, :, :, sIndex) =  ...
            - 1/2*(xTotalFluxStar(2:end, :, :) - xTotalFluxStar(1:end-1, :, :)) ./ dxC(2:end-1) ...
            - 1/2*(xTotalFluxOld(2:end, :, :) - xTotalFluxOld(1:end-1, :, :)) ./ dxC(2:end-1);
    else
        implicitTransportRHSTerms(:, :, :, sIndex) =  ...
            - (xTotalFluxStar(2:end, :, :) - xTotalFluxStar(1:end-1, :, :)) ./ dxC(2:end-1);
    end
        
%         + timeAlpha * (-(yTotalFluxOld(:, 2:end, :) - yTotalFluxOld(:, 1:end-1, :)) ./ reshape(dyC, [1, yN]) ...
%                        - (zTotalFluxOld(:, :, 2:end) - zTotalFluxOld(:, :, 1:end-1)) ./ reshape(dzC, [1, 1, zN])) ...
%         + timeBeta * (-(yTotalFluxOldOld(:, 2:end, :) - yTotalFluxOldOld(:, 1:end-1, :)) ./ reshape(dyC, [1, yN]) ...
%                        - (zTotalFluxOldOld(:, :, 2:end) - zTotalFluxOldOld(:, :, 1:end-1)) ./ reshape(dzC, [1, 1, zN])) ...
%         + timeGamma * (-(yTotalFluxOldOldOld(:, 2:end, :) - yTotalFluxOldOldOld(:, 1:end-1, :)) ./ reshape(dyC, [1, yN]) ...
%                        - (zTotalFluxOldOldOld(:, :, 2:end) - zTotalFluxOldOldOld(:, :, 1:end-1)) ./ reshape(dzC, [1, 1, zN])) ...
%         + timeDelta * (-(yTotalFluxOldOldOldOld(:, 2:end, :) - yTotalFluxOldOldOldOld(:, 1:end-1, :)) ./ reshape(dyC, [1, yN]) ...
%                        - (zTotalFluxOldOldOldOld(:, :, 2:end) - zTotalFluxOldOldOldOld(:, :, 1:end-1)) ./ reshape(dzC, [1, 1, zN]));
end

% Compute right hand side - Faradaic reactions
for fRxnIndex = 1:length(faradaicRxnInfo)
    if faradaicRxnInfo(fRxnIndex).reactants ~= 0
        individualFaradaicReactionTerm = 1/NA*faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
            * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
            * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - fields3DStar(2*nVar:nVar:end-nVar, 2:end-1, 2:end-1) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
            .* fields3DStar(nVar+faradaicRxnInfo(fRxnIndex).reactants:nVar:end-nVar, 2:end-1, 2:end-1)/faradaicRxnInfo(fRxnIndex).cRef;
        
%         individualFaradaicReactionTermDelta = 1/NA*faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
%             * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
%             * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - fields3DStar(2*nVar:nVar:end-nVar, 2:end-1, 2:end-1) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
%             .* fields3DStar(nVar+faradaicRxnInfo(fRxnIndex).reactants:nVar:end-nVar, 2:end-1, 2:end-1)/faradaicRxnInfo(fRxnIndex).cRef ...
%             .* (-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(-deltaFields(2*nVar:nVar:end-nVar, 2:end-1, 2:end-1)));

        faradaicRxnSourceTerms(:, :, :, faradaicRxnInfo(fRxnIndex).reactants) = ...
            faradaicRxnSourceTerms(:, :, :, faradaicRxnInfo(fRxnIndex).reactants) ...
            + individualFaradaicReactionTerm/e*faradaicRxnInfo(fRxnIndex).moleReactantsPerMoleElectron;% ...
            %+ individualFaradaicReactionTermDelta/e*faradaicRxnInfo(fRxnIndex).moleReactantsPerMoleElectron;

        if faradaicRxnInfo(fRxnIndex).products(1) ~= 0
            faradaicRxnSourceTerms(:, :, :, faradaicRxnInfo(fRxnIndex).products(1)) = ...
                faradaicRxnSourceTerms(:, :, :, faradaicRxnInfo(fRxnIndex).products(1)) ...
                - individualFaradaicReactionTerm/e*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(1); ...
                %- individualFaradaicReactionTermDelta/e*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(1);
        end

        if faradaicRxnInfo(fRxnIndex).products(2) ~= 0
            faradaicRxnSourceTerms(:, :, :, faradaicRxnInfo(fRxnIndex).products(2)) = ...
                faradaicRxnSourceTerms(:, :, :, faradaicRxnInfo(fRxnIndex).products(2)) ...
                - individualFaradaicReactionTerm/e*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(2); ...
                %- individualFaradaicReactionTermDelta/e*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(2);
        end
    else
        individualFaradaicReactionTerm = 1/NA*faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
            * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
            * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - fields3DStar(2*nVar:nVar:end-nVar, 2:end-1, 2:end-1) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential));
        
%         individualFaradaicReactionTermDelta = 1/NA*faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
%             * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
%             * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - fields3DStar(2*nVar:nVar:end-nVar, 2:end-1, 2:end-1) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
%             .* (-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(-deltaFields(2*nVar:nVar:end-nVar, 2:end-1, 2:end-1)));

        if faradaicRxnInfo(fRxnIndex).products(1) ~= 0
            faradaicRxnSourceTerms(:, :, :, faradaicRxnInfo(fRxnIndex).products(1)) = ...
                faradaicRxnSourceTerms(:, :, :, faradaicRxnInfo(fRxnIndex).products(1)) ...
                - individualFaradaicReactionTerm/e*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(1); ...
                %- individualFaradaicReactionTermDelta/e*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(1);
        end

        if faradaicRxnInfo(fRxnIndex).products(2) ~= 0
            faradaicRxnSourceTerms(:, :, :, faradaicRxnInfo(fRxnIndex).products(2)) = ...
                faradaicRxnSourceTerms(:, :, :, faradaicRxnInfo(fRxnIndex).products(2)) ...
                - individualFaradaicReactionTerm/e*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(2); ...
                %- individualFaradaicReactionTermDelta/e*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(2);
        end
    end
end

% Homogeneous reaction terms
rxns = constants.rxns{1};
for rxnIndex = 1:size(rxns, 1)
    if ((rxns(rxnIndex,1) ~= 0) && (rxns(rxnIndex,2) ~= 0))

        individualReactionTerm = -poro(1)*rxns(rxnIndex,5) ...
            * fields3DStar(nVar+rxns(rxnIndex,1):nVar:end-nVar, 2:end-1, 2:end-1) ...
            .* fields3DStar(nVar+rxns(rxnIndex,2):nVar:end-nVar, 2:end-1, 2:end-1);

        homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 1)) = homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 1)) + individualReactionTerm;

        homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 2)) = homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 2)) + individualReactionTerm;
        
        if rxns(rxnIndex,3) ~= 0
            homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 3)) = homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 3)) - individualReactionTerm;
        end

        if rxns(rxnIndex,4) ~= 0
            homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 4)) = homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 4)) - individualReactionTerm;
        end

    elseif ((rxns(rxnIndex,1) ~= 0) && (rxns(rxnIndex,2) == 0))
        individualReactionTerm = -poro(1)*rxns(rxnIndex,5) ...
            * fields3DStar(nVar+rxns(rxnIndex,1):nVar:end-nVar, 2:end-1, 2:end-1);
        
        homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 1)) = homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 1)) + individualReactionTerm;

        if rxns(rxnIndex,3) ~= 0
            homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 3)) = homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 3)) - individualReactionTerm;
        end

        if rxns(rxnIndex,4) ~= 0
            homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 4)) = homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 4)) - individualReactionTerm;
        end

    elseif ((rxns(rxnIndex,1) == 0) && (rxns(rxnIndex,2) ~= 0))
        individualReactionTerm = -poro(1)*rxns(rxnIndex,5) ...
            * fields3DStar(nVar+rxns(rxnIndex,2):nVar:end-nVar, 2:end-1, 2:end-1);
        
        homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 2)) = homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 2)) + individualReactionTerm;

        if rxns(rxnIndex,3) ~= 0
            homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 3)) = homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 3)) - individualReactionTerm;
        end

        if rxns(rxnIndex,4) ~= 0
            homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 4)) = homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 4)) - individualReactionTerm;
        end
    else
        % Add constant term to RHS
        individualReactionTerm = poro(1)*rxns(rxnIndex,5);
        if rxns(rxnIndex,3) ~= 0
            homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 3)) = homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 3)) + individualReactionTerm;
        end

        if rxns(rxnIndex,4) ~= 0
            homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 4)) = homogeneousRxnSourceTerms(:, :, :, rxns(rxnIndex, 4)) + individualReactionTerm;
        end
    end
end
    
end

%------------- END OF CODE --------------
