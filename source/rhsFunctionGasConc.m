function rhsVector = rhsFunctionGasConc(gasFields, gasVelocity, fields3D, dxC, dxF, dyC, dyF, dzC, dzF, nX, nY, nZ, gasLayerInfo, constants, uniqueSpecies, faradaicRxnInfo)

    nSpecies = length(uniqueSpecies);
    nGasSpecies = length(gasLayerInfo.gasSpeciesNames);
    
    rhsVector = zeros(size(gasFields));
    
    % One sided difference parameters
    deltaPlusZ = dzF(end-1) + dzF(end-2);
    gammaZ = 1/((deltaPlusZ)^2/dzF(end-1) - deltaPlusZ);
    betaZ = -deltaPlusZ^2/dzF(end-1)^2 * gammaZ;
    alphaZ = -betaZ - gammaZ;
    
    % advection component
    for sIndex = 1:nGasSpecies
        
        % Advection components
        speciesBCValue = gasLayerInfo.gasInitialMoleFraction(sIndex);
        concInterp = [speciesBCValue; 1/2*(gasFields(sIndex:nGasSpecies:end-nGasSpecies)+gasFields(nGasSpecies+sIndex:nGasSpecies:end))];
        
        [speciesFaradaicReactionRate, speciesDissolutionRate] = computeGasProductionRate(fields3D, dxC, dxF, dyC, dyF, dzC, dzF, nX, nY, nZ, gasLayerInfo, constants, uniqueSpecies, faradaicRxnInfo, sIndex);
        
        rhsVector(sIndex:nGasSpecies:end) = -1 ...
            * [1./dzC(1:end-1).*(concInterp(2:end).*gasVelocity(2:end-1) - concInterp(1:end-1).*gasVelocity(1:end-2)); ...
            (alphaZ * gasFields(end-nGasSpecies+sIndex) * 1/2*(gasVelocity(end) + gasVelocity(end-1)) ...
            + betaZ * gasFields(end-2*nGasSpecies+sIndex) * 1/2*(gasVelocity(end-1) + gasVelocity(end-2)) ...
            + gammaZ * gasFields(end-3*nGasSpecies+sIndex) * 1/2*(gasVelocity(end-2) + gasVelocity(end-3)))] ...
            + (speciesFaradaicReactionRate + speciesDissolutionRate) * (constants.NA * constants.kb * constants.T)/gasLayerInfo.gasPressure;
    end
    
end