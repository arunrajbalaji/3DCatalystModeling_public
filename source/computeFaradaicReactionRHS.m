function [rhsTerm] = computeFaradaicReactionRHS(xN, yN, zN, nSpecies, fields3DStar, faradaicRxnInfo, constants)
% Unpack constants and create convenience variables
phiElectrode = constants.phiElectrode;
NA = constants.NA;
kb = constants.kb;
T = constants.T;
e = constants.e;
nVar = nSpecies+1;

% Initialize variable for output, accumulating reactiont sources/sinks
rhsTerm = zeros(xN-2, yN, zN, nSpecies);

% Loop over individual reactions
% Compute right hand side - Faradaic reactions
for fRxnIndex = 1:length(faradaicRxnInfo)
    if faradaicRxnInfo(fRxnIndex).reactants ~= 0
        individualFaradaicReactionTerm = 1/NA*faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
            * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
            * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - fields3DStar(2*nVar:nVar:end-nVar, :, :) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
            .* fields3DStar(nVar+faradaicRxnInfo(fRxnIndex).reactants:nVar:end-nVar, :, :)/faradaicRxnInfo(fRxnIndex).cRef;
        
        rhsTerm(:, :, :, faradaicRxnInfo(fRxnIndex).reactants) = ...
            rhsTerm(:, :, :, faradaicRxnInfo(fRxnIndex).reactants) ...
            + individualFaradaicReactionTerm/e*faradaicRxnInfo(fRxnIndex).moleReactantsPerMoleElectron;

        if faradaicRxnInfo(fRxnIndex).products(1) ~= 0
            rhsTerm(:, :, :, faradaicRxnInfo(fRxnIndex).products(1)) = ...
                rhsTerm(:, :, :, faradaicRxnInfo(fRxnIndex).products(1)) ...
                - individualFaradaicReactionTerm/e*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(1);
        end

        if faradaicRxnInfo(fRxnIndex).products(2) ~= 0
            rhsTerm(:, :, :, faradaicRxnInfo(fRxnIndex).products(2)) = ...
                rhsTerm(:, :, :, faradaicRxnInfo(fRxnIndex).products(2)) ...
                - individualFaradaicReactionTerm/e*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(2);
        end
    else
        individualFaradaicReactionTerm = 1/NA*faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
            * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
            * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - fields3DStar(2*nVar:nVar:end-nVar, :, :) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential));

        if faradaicRxnInfo(fRxnIndex).products(1) ~= 0
            rhsTerm(:, :, :, faradaicRxnInfo(fRxnIndex).products(1)) = ...
                rhsTerm(:, :, :, faradaicRxnInfo(fRxnIndex).products(1)) ...
                - individualFaradaicReactionTerm/e*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(1);
        end

        if faradaicRxnInfo(fRxnIndex).products(2) ~= 0
            rhsTerm(:, :, :, faradaicRxnInfo(fRxnIndex).products(2)) = ...
                rhsTerm(:, :, :, faradaicRxnInfo(fRxnIndex).products(2)) ...
                - individualFaradaicReactionTerm/e*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(2);
        end
    end
end

end