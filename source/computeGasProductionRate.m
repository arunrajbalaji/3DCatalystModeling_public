function [speciesFaradaicReactionRate, speciesDissolutionRate] = computeGasProductionRate(fields3D, dxC, dxF, dyC, dyF, dzC, dzF, nX, nY, nZ, gasLayerInfo, constants, uniqueSpecies, faradaicRxnInfo, sIndex)

nSpecies = length(uniqueSpecies);

NA = constants.NA;
kb = constants.kb;
T = constants.T;
e = constants.e;
phiElectrode = constants.phiElectrode;
diff = constants.diff;
poro = constants.poro;
tort = constants.tort;

effectiveDiff = diff * poro(1) / tort(1);

% Coupling to electrolyte fluxes
speciesFaradaicReactionRate = zeros(length(dzC), 1);
speciesDissolutionRate = zeros(length(dzC), 1);

if gasLayerInfo.isTrackedInElectrolyte(sIndex)
    aqueousIndex = gasLayerInfo.isTrackedInElectrolyte(sIndex);
    
    % First, gas consumption due to dissolution at boundary, for
    % each slice in z
    for zIndex = 1:nZ
        % Compute flux at x=0 boundary
        % (neutral gas species ALWAYS have Dirichlet-type boundary
        % treatment, so flux must be computed using one-sided
        % scheme
        % NO electromigration for dissolved gaseous species

%         xFrontFlux = (- effectiveDiff(1,aqueousIndex) * (alphaXFront * fields3D(aqueousIndex,:,zIndex) ...
%             + betaXFront * fields3D((nSpecies+1)+aqueousIndex,:,zIndex) ...
%             + gammaXFront * fields3D(2*(nSpecies+1)+aqueousIndex,:,zIndex)))';

        xFrontFlux = (- effectiveDiff(1,aqueousIndex) * (2*fields3D((nSpecies+1)+aqueousIndex,:,zIndex) ...
            - 2*fields3D(aqueousIndex,:,zIndex)) / (dxF(1)))';
        
%         xBackFlux = (- effectiveDiff(1,aqueousIndex) * (alphaXBack * fields3D(end-(nSpecies+1)+aqueousIndex,:,zIndex) ...
%             + betaXBack * fields3D(end-2*(nSpecies+1)+aqueousIndex,:,zIndex) ...
%             + gammaXBack* fields3D(end-3*(nSpecies+1)+aqueousIndex,:,zIndex)))';

        xBackFlux = (- effectiveDiff(1,aqueousIndex) * (2*fields3D(end-(nSpecies+1)+aqueousIndex,:,zIndex) ...
            - 2*fields3D(end-2*(nSpecies+1)+aqueousIndex,:,zIndex)) / (dxF(end)))';
        
        % Integrate flux over 2 boundaries, to determine overall
        % consumption for cross sectional location in z
        integratedFluxInPlane = -sum(xFrontFlux .* dyC) + sum(xBackFlux .* dyC);
        
        % Convert consumption to correct units (per unit volume of gas)
        speciesDissolutionRate(zIndex) = (integratedFluxInPlane  * gasLayerInfo.nFingers) / gasLayerInfo.gasXYArea;
    end
    
elseif any(gasLayerInfo.isFaradaicProduct{sIndex})
    for fRxnIndex = 1:length(faradaicRxnInfo)
        if gasLayerInfo.isFaradaicProduct{sIndex}(fRxnIndex)
            for zIndex = 1:nZ
                if faradaicRxnInfo(fRxnIndex).reactants ~= 0
                    totalReactionRateInPlane = 0;
                    for xIndex = 2:nX-1
                        for yIndex = 1:nY
                            totalReactionRateInPlane = totalReactionRateInPlane ...
                                + dxC(xIndex)*dyC(yIndex) * 1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(gasLayerInfo.isFaradaicProduct{sIndex}(fRxnIndex))*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                                * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                                * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - fields3D(xIndex*(nSpecies+1), yIndex, zIndex) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                                * fields3D((xIndex-1)*(nSpecies+1) + faradaicRxnInfo(fRxnIndex).reactants, yIndex, zIndex)/faradaicRxnInfo(fRxnIndex).cRef;
                        end
                    end
                else
                    totalReactionRateInPlane = 0;
                    for xIndex = 2:nX-1
                        for yIndex = 1:nY
                            totalReactionRateInPlane = totalReactionRateInPlane ...
                                + dxC(xIndex)*dyC(yIndex) * 1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(gasLayerInfo.isFaradaicProduct{sIndex}(fRxnIndex))*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                                * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                                * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - fields3D(xIndex*(nSpecies+1), yIndex, zIndex) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential));
                        end
                    end
                end
                % ANSWER ACHIEVED FOR FIXED Z
                % add RHSTerm at THIS z to speciesFluxCellCenters
                speciesFaradaicReactionRate(zIndex) = (totalReactionRateInPlane * gasLayerInfo.nFingers) / gasLayerInfo.gasXYArea;
            end
        end
    end
end

end