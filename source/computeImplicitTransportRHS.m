function [rhsTerms] = computeImplicitTransportRHS(xN, yN, zN, dxC, dxF, nSpecies, fields3DStar, isFrontNoFluxCondition, isBackNoFluxCondition, constants, isNeutralGas)
% Unpack constants and create convenience variables
poro = constants.poro;
tort = constants.tort;
diff = constants.diff;
vale = constants.vale;
e = constants.e;
kb = constants.kb;
T = constants.T;
effectiveDiffVec = diff(1,:)' * poro(1) / tort(1);
nVar = nSpecies+1;

rhsTerms = zeros(xN-2, yN, zN, nSpecies);

for sIndex = 1:nSpecies
    % Compute right hand side - X Fluxes
    xDiff = -effectiveDiffVec(sIndex) * (fields3DStar(nVar+sIndex:nVar:end, :, :) - fields3DStar(sIndex:nVar:end-nVar, :, :)) ./ dxF;
    xElectro = -effectiveDiffVec(sIndex) * vale(1,sIndex)*e/(kb*T) * (fields3DStar(nVar+sIndex:nVar:end, :, :) + fields3DStar(sIndex:nVar:end-nVar, :, :))/2 ...
        .* (fields3DStar(2*nVar:nVar:end, :, :) - fields3DStar(nVar:nVar:end-nVar, :, :)) ./ dxF;
    
    xTotalFluxStar = xDiff + xElectro;

    if isFrontNoFluxCondition && ~isNeutralGas(sIndex)
        xTotalFluxStar(1, :, :) = 0;
    end

    if isBackNoFluxCondition && ~isNeutralGas(sIndex)
        xTotalFluxStar(end,: , :) = 0;
    end

    % Override for neutral gasses (correct Dirichlet BC)
    if isNeutralGas(sIndex)
        xTotalFluxStar(1, :, :) = 2 * effectiveDiffVec(sIndex)/dxF(1) * (fields3DStar(sIndex, :, :) ...
            - fields3DStar(nVar+sIndex, :, :));
        xTotalFluxStar(end, :, :) = 2 * effectiveDiffVec(sIndex)/dxF(end) *(-fields3DStar(end-nVar+sIndex, :, :) ...
            + fields3DStar(end-2*nVar+sIndex, :, :));
    end
    
    rhsTerms(:, :, :, sIndex) =  ...
        - (xTotalFluxStar(2:end, :, :) - xTotalFluxStar(1:end-1, :, :)) ./ dxC(2:end-1);
end

end