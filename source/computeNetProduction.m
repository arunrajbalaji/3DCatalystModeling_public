function netProductionRateStrip = computeNetProduction(dxC, dyC, dzC, nSpecies, xStripMiddle, constants, nEq, faradaicRxnInfo)

poro = constants.poro;
e = constants.e;
kb = constants.kb;
T = constants.T;
NA = constants.NA;
wienBeta = constants.wienBeta;
phiElectrode = constants.phiElectrode;

maxWienFactor = 0;

% Might need this later, implemented here in a very ad-hoc way. Leaving for
% reference, definitely improve if necessary.
% deltaStern = 1e-10;
% exchCurrent = 0.0;
% convFactor = 1/(1e3 / 1e4 * 6.022e23 * 1.602e-19);
% faradaicCoeff = deltaStern*e/kb/T;

xN = length(dxC);

b = zeros(nEq,1);
netProductionRateStrip = zeros(nSpecies,1);

% Add reaction terms as necessary
% Loop over grid points
for xLocIndex = 2:xN-1
    % Add reaction terms for each cell
    rxns = constants.rxns{xLocIndex};
    
    for rxnIndex = 1:size(rxns, 1)
        % Reactants side first
        if ((rxns(rxnIndex,1) ~= 0) && (rxns(rxnIndex,2) ~= 0))
            
            b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,1)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,1))...
                - poro(xLocIndex)*rxns(rxnIndex,5)...
                * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))...
                * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
            
            b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,2)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,2))...
                - poro(xLocIndex)*rxns(rxnIndex,5)...
                * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))...
                * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
            
            if rxns(rxnIndex,3) ~= 0
                
                b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3))...
                    + poro(xLocIndex)*rxns(rxnIndex,5)...
                    * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))...
                    * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
            end
            
            if rxns(rxnIndex,4) ~= 0
                
                b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4))...
                    + poro(xLocIndex)*rxns(rxnIndex,5)...
                    * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))...
                    * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
            end
            
        elseif ((rxns(rxnIndex,1) ~= 0) && (rxns(rxnIndex,2) == 0))
            dPhidXStar = ((xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1) + (nSpecies+1)) + xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)))/2 ...
                - (xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)) + xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1)))/2)/dxC(xLocIndex);
            wienExponential = exp(rxns(rxnIndex,6)*wienBeta*abs(dPhidXStar));
            
            if wienExponential > maxWienFactor && rxns(rxnIndex,6)~=0
                maxWienFactor = wienExponential;
            end
            
            b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,1)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,1))...
                - poro(xLocIndex)*rxns(rxnIndex,5) * wienExponential * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
            
            if rxns(rxnIndex,3) ~= 0
                
                b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3))...
                    + poro(xLocIndex)*rxns(rxnIndex,5) * wienExponential * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
                
            end
            
            if rxns(rxnIndex,4) ~= 0
                
                b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4))...
                    + poro(xLocIndex)*rxns(rxnIndex,5) * wienExponential * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
                
            end
            
        elseif ((rxns(rxnIndex,1) == 0) && (rxns(rxnIndex,2) ~= 0))
            wienExponential = exp(rxns(rxnIndex,6)*wienBeta*abs(((xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1) + (nSpecies+1)) + xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)))/2 ...
                - (xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1) + (nSpecies+1)) + xStripMiddle((xLocIndex-1 - 1)*(nSpecies+1) + (nSpecies+1)))/2)/dxC(xLocIndex)));
            
            b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,2)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,2))...
                - poro(xLocIndex)*rxns(rxnIndex,5) * wienExponential * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
            
            if rxns(rxnIndex,3) ~= 0
                
                b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3))...
                    + poro(xLocIndex)*rxns(rxnIndex,5) * wienExponential * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
            end
            
            if rxns(rxnIndex,4) ~= 0
                
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
        end
        
        if faradaicRxnInfo(fRxnIndex).products(1) ~= 0
            if faradaicRxnInfo(fRxnIndex).reactants ~= 0
                b((xLocIndex-2)*nSpecies + faradaicRxnInfo(fRxnIndex).products(1)) = b((xLocIndex-2)*nSpecies + faradaicRxnInfo(fRxnIndex).products(1)) ...
                    +  1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(1)*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                    * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                    * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                    * xStripMiddle((xLocIndex-1)*(nSpecies+1) + faradaicRxnInfo(fRxnIndex).reactants)/faradaicRxnInfo(fRxnIndex).cRef;
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
            else
                b((xLocIndex-2)*nSpecies + faradaicRxnInfo(fRxnIndex).products(2)) = b((xLocIndex-2)*nSpecies + faradaicRxnInfo(fRxnIndex).products(2)) ...
                    + 1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(2)*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                    * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                    * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential));
            end
        end
    end
end

for sIndex = 1:nSpecies
    netProductionRateStrip(sIndex) = dyC * dzC * sum(dxC .* b(sIndex:(nSpecies+1):end));
end

end