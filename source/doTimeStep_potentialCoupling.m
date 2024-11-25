function [deltaPhi, residualNorm] = doTimeStep_potentialCoupling(dxC, dxF, dyC, dyF, dzC, dzF, ...
    nSpecies, fields3DStar, isFrontNoFluxCondition, isBackNoFluxCondition, isLeftNoFluxCondition, isRightNoFluxCondition, isBottomNoFluxCondition, isTopNoFluxCondition, ...
    constants, yLeftBCVec, yRightBCVec, zBottomBCVec, zTopBCVec, faradaicRxnInfo, isNeutralGas, yRightElectrolyteOverride, waterRxnIndices, ohIndex, hIndex)
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

RHSMatrix = zeros(xN-2, yN, zN);
stericSum = zeros(xN, yN+2, zN+2);
RHSFaradaicReactions = zeros(xN-2, yN, zN);
dPhiCenterDiagFaradaicComponent = zeros(xN-2, yN, zN);

for xIndex = 1:xN
    stericSum(xIndex, :, :) = sum(fields3DStar((xIndex-1)*nVar+1:(xIndex-1)*nVar+nSpecies, :, :) .* stericOnOffVec, 1);
end

for sIndex = 1:nSpecies
    % Compute right hand side - X Fluxes
    xDiff = -effectiveDiffVec(sIndex) * (fields3DStar(nVar+sIndex:nVar:end, 2:end-1, 2:end-1) - fields3DStar(sIndex:nVar:end-nVar, 2:end-1, 2:end-1)) ./ dxF;
    xElectro = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * (fields3DStar(nVar+sIndex:nVar:end, 2:end-1, 2:end-1) + fields3DStar(sIndex:nVar:end-nVar, 2:end-1, 2:end-1))/2 ...
        .* (fields3DStar(2*nVar:nVar:end, 2:end-1, 2:end-1) - fields3DStar(nVar:nVar:end-nVar, 2:end-1, 2:end-1)) ./ dxF;
    xSteric = -effectiveDiffVec(sIndex) * (fields3DStar(nVar+sIndex:nVar:end, 2:end-1, 2:end-1) + fields3DStar(sIndex:nVar:end-nVar, 2:end-1, 2:end-1))/2 ...
        .* (log(1 - stericACubed * stericSum(2:end, 2:end-1, 2:end-1)) - log(1 - stericACubed * stericSum(1:end-1, 2:end-1, 2:end-1))) ./ dxF;

    xTotalFlux = xDiff + xElectro + xSteric;

    if isFrontNoFluxCondition && ~isNeutralGas(sIndex)
        xTotalFlux(1, :, :) = 0;
    end

    if isBackNoFluxCondition && ~isNeutralGas(sIndex)
        xTotalFlux(end,: , :) = 0;
    end

    if isFrontNoFluxCondition && ~isNeutralGas(sIndex)
        xTotalFluxDelta(1, :, :) = 0;
    end

    if isBackNoFluxCondition && ~isNeutralGas(sIndex)
        xTotalFluxDelta(end,: , :) = 0;
    end


    if yRightElectrolyteOverride
        if (sIndex == hIndex)   % HYDROGEN
            % Find OH- concentration at ghost cell center
            ghostCellOH = -dyC(end)/2/dyF(end-1)*fields3DStar(ohIndex:nVar:end, end-2, :) + (dyC(end)/2 + dyF(end-1))/(dyF(end-1)) * fields3DStar(ohIndex:nVar:end, end-1, :);
            % Equilibrate water recombination
            % Set Dirichlet value
            K_eq = constants.rxns{1}(waterRxnIndices(1),5)/constants.rxns{1}(waterRxnIndices(2),5);
            fields3DStar(hIndex:nVar:end, end, :) = K_eq ./ ghostCellOH;
        end
    end

    % Compute right hand side - Y Fluxes
    yDiff = -effectiveDiffVec(sIndex) * (fields3DStar(nVar+sIndex:nVar:end-nVar, 2:end, 2:end-1) - fields3DStar(nVar+sIndex:nVar:end-nVar, 1:end-1, 2:end-1)) ./ reshape(dyF, [1, yN+1]);
    yElectro = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * (fields3DStar(nVar+sIndex:nVar:end-nVar, 2:end, 2:end-1) + fields3DStar(nVar+sIndex:nVar:end-nVar, 1:end-1, 2:end-1))/2 ...
        .* (fields3DStar(2*nVar:nVar:end-nVar, 2:end, 2:end-1) - fields3DStar(2*nVar:nVar:end-nVar, 1:end-1, 2:end-1)) ./ reshape(dyF, [1, yN+1]);
    % Correct Dirichlet treatment for yRight boundary
    yDiff(:, end, :) = -effectiveDiffVec(sIndex) * 2 * (fields3DStar(nVar+sIndex:nVar:end-nVar, end, 2:end-1) - fields3DStar(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1)) ./ dyF(end);
    yElectro(:, end, :) = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * fields3DStar(nVar+sIndex:nVar:end-nVar, end, 2:end-1) ...
        * 2 .* (fields3DStar(2*nVar:nVar:end-nVar, end, 2:end-1) - fields3DStar(2*nVar:nVar:end-nVar, end-1, 2:end-1)) ./ dyF(end);
    ySteric = -effectiveDiffVec(sIndex) * (fields3DStar(nVar+sIndex:nVar:end-nVar, 2:end, 2:end-1) + fields3DStar(nVar+sIndex:nVar:end-nVar, 1:end-1, 2:end-1))/2 ...
        .* (log(1 - stericACubed * stericSum(2:end-1, 2:end, 2:end-1)) - log(1 - stericACubed * stericSum(2:end-1, 1:end-1, 2:end-1))) ./ reshape(dyF, [1, yN+1]);

    if yRightElectrolyteOverride
        if vale(1, sIndex) < 0
            delta_1_2 = dyC(end)/2;
            delta_3_2 = dyC(end)/2 + dyF(end-1);
            delta_5_2 = dyC(end)/2 + dyF(end-1) + dyF(end-2);
            beta = (-1/(delta_5_2 - delta_1_2)*(delta_1_2^2 - delta_5_2^2)) / (-(delta_5_2 - delta_3_2)/(delta_5_2 - delta_1_2) * (delta_1_2^2 - delta_5_2^2) + delta_3_2^2 - delta_5_2^2);
            alpha = 1/(delta_5_2 - delta_1_2) - beta * (delta_5_2 - delta_3_2)/(delta_5_2 - delta_1_2);
            gamma = -alpha - beta;
            
            yDiff(:, end, :) = - effectiveDiffVec(sIndex) * (alpha * fields3DStar(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1) + beta * fields3DStar(nVar+sIndex:nVar:end-nVar, end-2, 2:end-1) + gamma * fields3DStar(nVar+sIndex:nVar:end-nVar, end-3, 2:end-1));
            yElectro(:, end, :) = -vale(1,sIndex)*e/(kb*T)*effectiveDiffVec(sIndex) * (-dyC(end)/2/dyF(end-1) * fields3DStar(nVar+sIndex:nVar:end-nVar, end-2, 2:end-1) + (dyC(end)/2 + dyF(end-1))/(dyF(end-1)) * fields3DStar(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1)) ...
                * 2 .* (fields3DStar(2*nVar:nVar:end-nVar, end, 2:end-1) - fields3DStar(2*nVar:nVar:end-nVar, end-1, 2:end-1))/dyF(end);
            effectiveBCConc = -dyC(end)/2/dyF(end-1) * fields3DStar(:, end-2, :) + (dyC(end)/2 + dyF(end-1))/(dyF(end-1)) * fields3DStar(:, end-1, :);
            
            fields3DStar(sIndex:nVar:end, end, :) = effectiveBCConc(sIndex:nVar:end, 1, :);
            stericSum(:, end, :) = 0;
            for xIndex = 1:xN
                stericSum(:, end, :) = stericSum(:, end, :) + sum(effectiveBCConc((xIndex-1)*nVar+1:(xIndex-1)*nVar+nSpecies, 1, :) .* stericOnOffVec, 1);
            end
            
            ySteric(:, end, :) = -effectiveDiffVec(sIndex) * (fields3DStar(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1) + effectiveBCConc(nVar+sIndex:nVar:end-nVar, 1, 2:end-1))/2 ...
                .* (log(1 - stericACubed * stericSum(2:end-1, end, 2:end-1)) - log(1 - stericACubed * stericSum(2:end-1, end-1, 2:end-1))) ./ dyF(end);
            % Commenting this out ensures that we are using simple
            % Dirichlet for neutrals and cations.
%         elseif vale(1, sIndex) == 0
%             yDiff(:, end, :) = -poro(1)/tort(1)*(diff(1,1)+diff(1,4))/2 / bulkElectrolyteL * (bulkElectrolyteC(sIndex) - (-1/2 * fields3DStar(nVar+sIndex:nVar:end-nVar, end-2, 2:end-1) + 3/2 * fields3DStar(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1)));
%             yElectro(:, end, :) = 0;
%             ySteric(:, end, :) = 0;
% 
%         elseif vale(1, sIndex) > 0
%             effectiveBCConc = -1 * fields3DStar(:, end-2, :) + 2 * fields3DStar(:, end-1, :);
%             stericSum(:, end, :) = 0;
%             for xIndex = 1:xN
%                 stericSum(:, end, :) = stericSum(:, end, :) + sum(effectiveBCConc((xIndex-1)*nVar+1:(xIndex-1)*nVar+nSpecies, 1, :) .* stericOnOffVec, 1);
%             end
%             ySteric(:, end, :) = -effectiveDiffVec(sIndex) * (fields3DStar(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1) + effectiveBCConc(nVar+sIndex:nVar:end-nVar, 1, 2:end-1))/2 ...
%                 .* (log(1 - stericACubed * stericSum(2:end-1, end, 2:end-1)) - log(1 - stericACubed * stericSum(2:end-1, end-1, 2:end-1))) ./ dyF(end);
%             yElectro(:, end, :) = -vale(1,sIndex)*e/(kb*T)*effectiveDiffVec(sIndex) * (fields3DStar(nVar+sIndex:nVar:end-nVar, end, 2:end-1) + fields3DStar(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1))/2 ...
%                 * 2 .* (fields3DStar(2*nVar:nVar:end-nVar, end, 2:end-1) - fields3DStar(2*nVar:nVar:end-nVar, end-1, 2:end-1))/dyF(end);
        end
    end

    yTotalFlux = yDiff + yElectro + ySteric;

    if isLeftNoFluxCondition && ~isNeutralGas(sIndex)
        yTotalFlux(:, 1, :) = 0;
    end

    if isRightNoFluxCondition && ~isNeutralGas(sIndex) && ~yRightElectrolyteOverride
        yTotalFlux(:, end, :) = 0;
    end

    % Override for neutral gasses (correct Dirichlet BC) OFF FOR NOW
    if isNeutralGas(sIndex)
        %         yTotalFluxOld(:, 1, :) = 2 * effectiveDiffVec(sIndex)/dyF(1) * fields3DOld(nVar+sIndex:nVar:end-nVar, 1, 2:end-1) ...
        %             - 2 * effectiveDiffVec(sIndex)/dyF(1) * fields3DOld(nVar+sIndex:nVar:end-nVar, 2, 2:end-1);
        yTotalFlux(:, 1, :) = 0;
    end

    % Compute right hand side - Z Fluxes
    zDiff = -effectiveDiffVec(sIndex) * (fields3DStar(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end) - fields3DStar(nVar+sIndex:nVar:end-nVar, 2:end-1, 1:end-1)) ./ reshape(dzF, [1, 1, zN+1]);
    zElectro = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * (fields3DStar(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end) + fields3DStar(nVar+sIndex:nVar:end-nVar, 2:end-1, 1:end-1))/2 ...
        .* (fields3DStar(2*nVar:nVar:end-nVar, 2:end-1, 2:end) - fields3DStar(2*nVar:nVar:end-nVar, 2:end-1, 1:end-1)) ./ reshape(dzF, [1, 1, zN+1]);
    zSteric = -effectiveDiffVec(sIndex) * (fields3DStar(nVar+sIndex:nVar:end-nVar, 2:end-1, 2:end) + fields3DStar(nVar+sIndex:nVar:end-nVar, 2:end-1, 1:end-1))/2 ...
        .* (log(1 - stericACubed * stericSum(2:end-1, 2:end-1, 2:end)) - log(1 - stericACubed * stericSum(2:end-1, 2:end-1, 1:end-1))) ./ reshape(dzF, [1, 1, zN+1]);


    zTotalFlux = zDiff + zElectro + zSteric;

    if isBottomNoFluxCondition
        zTotalFlux(:, :, 1) = 0;
    end

    if isTopNoFluxCondition
        zTotalFlux(:, :, end) = 0;
    end


    % Flux has a negative in the definition
    % Flux appears on RHS with negative in front (move over)
    RHSMatrix = RHSMatrix ...
        - e * vale(1, sIndex) * (xTotalFlux(2:end, :, :) - xTotalFlux(1:end-1, :, :)) ./ dxC(2:end-1) ...
        - e * vale(1, sIndex) * (yTotalFlux(:, 2:end, :) - yTotalFlux(:, 1:end-1, :)) ./ reshape(dyC, [1, yN]) ...
        - e * vale(1, sIndex) * (zTotalFlux(:, :, 2:end) - zTotalFlux(:, :, 1:end-1)) ./ reshape(dzC, [1, 1, zN]);

end

% Compute right hand side - Faradaic reactions
for fRxnIndex = 1:length(faradaicRxnInfo)
    if faradaicRxnInfo(fRxnIndex).reactants ~= 0
        individualFaradaicReactionTerm = 1/NA*faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
            * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
            * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - fields3DStar(2*nVar:nVar:end-nVar, 2:end-1, 2:end-1) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
            .* fields3DStar(nVar+faradaicRxnInfo(fRxnIndex).reactants:nVar:end-nVar, 2:end-1, 2:end-1)/faradaicRxnInfo(fRxnIndex).cRef;

        RHSFaradaicReactions = RHSFaradaicReactions + individualFaradaicReactionTerm;% + individualFaradaicReactionTermDelta;
    else
        individualFaradaicReactionTerm = 1/NA*faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
            * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
            * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - fields3DStar(2*nVar:nVar:end-nVar, 2:end-1, 2:end-1) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential));
        RHSFaradaicReactions = RHSFaradaicReactions + individualFaradaicReactionTerm;
    end

    dPhiCenterDiagFaradaicComponent = dPhiCenterDiagFaradaicComponent ...
        - individualFaradaicReactionTerm * faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T);
end
RHSMatrix = RHSMatrix + RHSFaradaicReactions;


electroSum = zeros(xN, yN+2, zN+2);
for xIndex = 1:xN
    electroSum(xIndex, :, :) = sum((vale(1,:).^2)' .* effectiveDiffVec .* fields3DStar((xIndex-1)*nVar+1:(xIndex-1)*nVar+nSpecies, :, :), 1);
end

% Compute matrix entries
dPhiCenterDiagX = -e^2/(kb*T)./dxC(2:end-1) ...
    .* [-(electroSum(2, 2:end-1, 2:end-1) + electroSum(3, 2:end-1, 2:end-1))/2./dxF(2); ...
    (-(electroSum(3:end-2, 2:end-1, 2:end-1) + electroSum(4:end-1, 2:end-1, 2:end-1))/2./dxF(3:end-1) ...
    - (electroSum(3:end-2, 2:end-1, 2:end-1) + electroSum(2:end-3, 2:end-1, 2:end-1))/2./dxF(2:end-2)); ...
    -(electroSum(end-1, 2:end-1, 2:end-1) + electroSum(end-2, 2:end-1, 2:end-1))/2./dxF(end-1)];

dPhiCenterDiagY = -e^2/(kb*T)./reshape(dyC(1:end), [1, yN]) ...
    .* cat(2, -(electroSum(2:end-1, 2, 2:end-1) + electroSum(2:end-1, 3, 2:end-1))/2/dyF(2), ...
    (-(electroSum(2:end-1, 3:end-2, 2:end-1) + electroSum(2:end-1, 2:end-3, 2:end-1))/2./reshape(dyF(2:end-2), [1, length(dyF)-3]) ...
    - (electroSum(2:end-1, 3:end-2, 2:end-1) + electroSum(2:end-1, 4:end-1, 2:end-1))/2./reshape(dyF(3:end-1), [1, length(dyF)-3])), ...
    -(electroSum(2:end-1, end-1, 2:end-1) + electroSum(2:end-1, end-2, 2:end-1))/2./dyF(end-1) ...
    - 2*electroSum(2:end-1, end, 2:end-1)./dyF(end));

dPhiCenterDiagZ = -e^2/(kb*T)./reshape(dzC(1:end), [1, 1, zN]) ...
    .* cat(3, -(electroSum(2:end-1, 2:end-1, 2) + electroSum(2:end-1, 2:end-1, 3))/2/dzF(2), ...
    (-(electroSum(2:end-1, 2:end-1, 3:end-2) + electroSum(2:end-1, 2:end-1, 2:end-3))/2./reshape(dzF(2:end-2), [1, 1, length(dzF)-3]) ...
    - (electroSum(2:end-1, 2:end-1, 3:end-2) + electroSum(2:end-1, 2:end-1, 4:end-1))/2./reshape(dzF(3:end-1), [1, 1, length(dzF)-3])), ...
    -(electroSum(2:end-1, 2:end-1, end-1) + electroSum(2:end-1, 2:end-1, end-2))/2./dzF(end-1));

dPhiCenterDiag = dPhiCenterDiagX  + dPhiCenterDiagFaradaicComponent + dPhiCenterDiagY + dPhiCenterDiagZ;% + dPhiCenterDiagXDelta;

dPhiFrontDiag = cat(1, zeros(1, yN, zN), -e^2/(kb*T)./dxC(3:end-1) .* (electroSum(3:end-1, 2:end-1, 2:end-1) + electroSum(2:end-2, 2:end-1, 2:end-1))/2./dxF(2:end-1));
dPhiBackDiag = cat(1, -e^2/(kb*T)./dxC(2:end-2) .* (electroSum(2:end-2, 2:end-1, 2:end-1) + electroSum(3:end-1, 2:end-1, 2:end-1))/2./dxF(2:end-1), zeros(1, yN, zN));

dPhiLeftDiag = cat(2, zeros(xN-2, 1, zN), -e^2/(kb*T)./reshape(dyC(2:end), [1, yN-1]) .* (electroSum(2:end-1, 3:end-1, 2:end-1) + electroSum(2:end-1, 2:end-2, 2:end-1))/2./reshape(dyF(2:end-1), [1, length(dyF)-2]));
dPhiRightDiag = cat(2, -e^2/(kb*T)./reshape(dyC(1:end-1), [1, yN-1]) .* (electroSum(2:end-1, 2:end-2, 2:end-1) + electroSum(2:end-1, 3:end-1, 2:end-1))/2./reshape(dyF(2:end-1), [1, length(dyF)-2]), zeros(xN-2, 1, zN));

dPhiBottomDiag = cat(3, zeros(xN-2, yN, 1), -e^2/(kb*T)./reshape(dzC(2:end), [1, 1, zN-1]) .* (electroSum(2:end-1, 2:end-1, 3:end-1) + electroSum(2:end-1, 2:end-1, 2:end-2))/2./reshape(dzF(2:end-1), [1, 1, length(dzF)-2]));
dPhiTopDiag = cat(3, -e^2/(kb*T)./reshape(dzC(1:end-1), [1, 1, zN-1]) .* (electroSum(2:end-1, 2:end-1, 2:end-2) + electroSum(2:end-1, 2:end-1, 3:end-1))/2./reshape(dzF(2:end-1), [1, 1, length(dzF)-2]), zeros(xN-2, yN, 1));

nMax = (xN-2)*yN*zN;

LHS = spdiags([...
    [dPhiBottomDiag((xN-2)*yN+1:nMax)'; zeros((xN-2)*yN, 1)], ...
    [dPhiLeftDiag(xN-1:nMax)'; zeros(xN-2, 1)], ...
    [(dPhiFrontDiag(2:nMax))'; 0], ...
    dPhiCenterDiag(1:nMax)', ...
    [0; (dPhiBackDiag(1:nMax-1))'], ...
    [zeros((xN-2),1); dPhiRightDiag(1:nMax-(xN-2))'], ...
    [zeros((xN-2)*yN, 1); dPhiTopDiag(1:nMax-(xN-2)*yN)']], ...
    ...
    [-(xN-2)*yN -(xN-2) -1 0 1 (xN-2) (xN-2)*yN], nMax, nMax);

% options = struct("type","ilutp","droptol",1e-6);    
% [L,U] = ilu(LHS, options);
% [deltaPhi, ~] = bicgstab(LHS, RHSMatrix(1:nMax)', 1e-6, 100, L, U);

warning('off', 'MATLAB:nearlySingularMatrix')
warning('off', 'MATLAB:singularMatrix')
deltaPhi = LHS \ RHSMatrix(1:nMax)';

deltaPhi = reshape(deltaPhi, [xN-2, yN, zN]);
deltaPhi = [zeros(1, yN, zN); deltaPhi; zeros(1, yN, zN)];

residualNorm = sqrt(sum(RHSMatrix.^2, [1 2 3]));

end

%------------- END OF CODE --------------
