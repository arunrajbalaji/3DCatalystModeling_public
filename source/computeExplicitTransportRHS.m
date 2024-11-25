function [rhsTerms] = computeExplicitTransportRHS(dyC, dyF, dzC, dzF, xN, yN, zN, fields3DStar, constants, ...
    isLeftNoFluxCondition, isRightNoFluxCondition, isBottomNoFluxCondition, isTopNoFluxCondition, ...
    yLeftBCVec, yRightBCVec, zBottomBCVec, zTopBCVec, yRightElectrolyteOverride, isNeutralGas, nSpecies, waterRxnIndices, ohIndex, hIndex)

poro = constants.poro;
tort = constants.tort;
diff = constants.diff;
vale = constants.vale;
e = constants.e;
kb = constants.kb;
T = constants.T;

nVar = nSpecies+1;

fields3DStar = cat(2, yLeftBCVec, fields3DStar, yRightBCVec);
fields3DStar = cat(3, cat(2, zBottomBCVec(:, 1, 1), zBottomBCVec, zBottomBCVec(:, end, 1)), ...
    fields3DStar, ...
    cat(2, zTopBCVec(:, 1, end), zTopBCVec, zTopBCVec(:, end, end)));

effectiveDiffVec = diff(1,:)' * poro(1) / tort(1);

rhsTerms = zeros(xN-2, yN, zN, nSpecies);

for sIndex = 1:nVar-1
    
    if yRightElectrolyteOverride
        if (sIndex == hIndex)   % HYDROGEN
            K_eq = constants.rxns{1}(waterRxnIndices(1),5)/constants.rxns{1}(waterRxnIndices(2),5);
            % Find OH- concentration at ghost cell center
            ghostCellOH = -dyC(end)/2/dyF(end-1)*fields3DStar(ohIndex:nVar:end, end-2, :) + (dyC(end)/2 + dyF(end-1))/(dyF(end-1)) * fields3DStar(ohIndex:nVar:end, end-1, :);
            % Equilibrate water recombinationm, Set Dirichlet value
            fields3DStar(hIndex:nVar:end, end, :) = K_eq ./ ghostCellOH;
        end
    end
    
    % Compute right hand side - Y Fluxes
    yDiff = -effectiveDiffVec(sIndex) * (fields3DStar(nVar+sIndex:nVar:end-nVar, 2:end, 2:end-1) - fields3DStar(nVar+sIndex:nVar:end-nVar, 1:end-1, 2:end-1)) ./ reshape(dyF, [1, yN+1]);
    yElectro = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * (fields3DStar(nVar+sIndex:nVar:end-nVar, 2:end, 2:end-1) + fields3DStar(nVar+sIndex:nVar:end-nVar, 1:end-1, 2:end-1))/2 ...
        .* (fields3DStar(2*nVar:nVar:end-nVar, 2:end, 2:end-1) - (fields3DStar(2*nVar:nVar:end-nVar, 1:end-1, 2:end-1))) ./ reshape(dyF, [1, yN+1]);
    % Correct Dirichlet treatment for yRight boundary
    yDiff(:, end, :) = -effectiveDiffVec(sIndex) * 2 * (fields3DStar(nVar+sIndex:nVar:end-nVar, end, 2:end-1) - fields3DStar(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1)) ./ dyF(end);
    yElectro(:, end, :) = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * fields3DStar(nVar+sIndex:nVar:end-nVar, end, 2:end-1) ...
        * 2 .* (fields3DStar(2*nVar:nVar:end-nVar, end, 2:end-1) - fields3DStar(2*nVar:nVar:end-nVar, end-1, 2:end-1)) ./ dyF(end);
    
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
                * 2 .* (fields3DStar(2*nVar:nVar:end-nVar, end, 2:end-1) - (fields3DStar(2*nVar:nVar:end-nVar, end-1, 2:end-1) ))/dyF(end);
            effectiveBCConc = -dyC(end)/2/dyF(end-1) * fields3DStar(:, end-2, :) + (dyC(end)/2 + dyF(end-1))/(dyF(end-1)) * fields3DStar(:, end-1, :);
            fields3DStar(sIndex:nVar:end, end, :) = effectiveBCConc(sIndex:nVar:end, 1, :);
            % Commenting this out ensures that we are using simple
            % Dirichlet for neutrals and cations.
%         elseif vale(1, sIndex) == 0
%             yDiff(:, end, :) = -poro(1)/tort(1)*(diff(1,1)+diff(1,4))/2 / bulkElectrolyteL * (bulkElectrolyteC(sIndex) - (-1/2 * fields3DStar(nVar+sIndex:nVar:end-nVar, end-2, 2:end-1) + 3/2 * fields3DStar(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1)));
%             yElectro(:, end, :) = 0;
%         elseif vale(1, sIndex) > 0
%             yElectro(:, end, :) = -vale(1,sIndex)*e/(kb*T)*effectiveDiffVec(sIndex) * (fields3DStar(nVar+sIndex:nVar:end-nVar, end, 2:end-1) + fields3DStar(nVar+sIndex:nVar:end-nVar, end-1, 2:end-1))/2 ...
%                 * 2 .* (fields3DStar(2*nVar:nVar:end-nVar, end, 2:end-1) - (fields3DStar(2*nVar:nVar:end-nVar, end-1, 2:end-1)))/dyF(end);
        end
    end
    
    yTotalFlux = yDiff + yElectro;
    
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
        .* (fields3DStar(2*nVar:nVar:end-nVar, 2:end-1, 2:end) - (fields3DStar(2*nVar:nVar:end-nVar, 2:end-1, 1:end-1))) ./ reshape(dzF, [1, 1, zN+1]);
    
    
    zTotalFlux = zDiff + zElectro;
    
    if isBottomNoFluxCondition
        zTotalFlux(:, :, 1) = 0;
    end
    
    if isTopNoFluxCondition
        zTotalFlux(:, :, end) = 0;
    end
    
    rhsTerms(:, :, :, sIndex) =  ...
        - (yTotalFlux(:, 2:end, :) - yTotalFlux(:, 1:end-1, :)) ./ reshape(dyC, [1, yN]) ...
        - (zTotalFlux(:, :, 2:end) - zTotalFlux(:, :, 1:end-1)) ./ reshape(dzC, [1, 1, zN]);
end

end