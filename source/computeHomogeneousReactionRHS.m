function [rhsTerm] = computeHomogeneousReactionRHS(xN, yN, zN, nSpecies, fields3DStar, constants)
% Unpack constants and create convenience variables
poro = constants.poro;
rxns = constants.rxns{1};
nVar = nSpecies+1;

% Initialize variable for output, accumulating reactiont sources/sinks
rhsTerm = zeros(xN-2, yN, zN, nSpecies);

% Loop over individual reactions
for rxnIndex = 1:size(rxns, 1)
    % Two reactants are present
    if ((rxns(rxnIndex,1) ~= 0) && (rxns(rxnIndex,2) ~= 0))

        individualReactionTerm = -poro(1)*rxns(rxnIndex,5) ...
            * fields3DStar(nVar+rxns(rxnIndex,1):nVar:end-nVar, :, :) ...
            .* fields3DStar(nVar+rxns(rxnIndex,2):nVar:end-nVar, :, :);

        rhsTerm(:, :, :, rxns(rxnIndex, 1)) = rhsTerm(:, :, :, rxns(rxnIndex, 1)) + individualReactionTerm;

        rhsTerm(:, :, :, rxns(rxnIndex, 2)) = rhsTerm(:, :, :, rxns(rxnIndex, 2)) + individualReactionTerm;
        
        % First product (with two reactants)
        if rxns(rxnIndex,3) ~= 0
            rhsTerm(:, :, :, rxns(rxnIndex, 3)) = rhsTerm(:, :, :, rxns(rxnIndex, 3)) - individualReactionTerm;
        end

        % Second product (with two reactants)
        if rxns(rxnIndex,4) ~= 0
            rhsTerm(:, :, :, rxns(rxnIndex, 4)) = rhsTerm(:, :, :, rxns(rxnIndex, 4)) - individualReactionTerm;
        end

    % Only first reactant is present
    elseif ((rxns(rxnIndex,1) ~= 0) && (rxns(rxnIndex,2) == 0))
        individualReactionTerm = -poro(1)*rxns(rxnIndex,5) ...
            * fields3DStar(nVar+rxns(rxnIndex,1):nVar:end-nVar, :, :);
        
        rhsTerm(:, :, :, rxns(rxnIndex, 1)) = rhsTerm(:, :, :, rxns(rxnIndex, 1)) + individualReactionTerm;

        % First product (with first reactant only)
        if rxns(rxnIndex,3) ~= 0
            rhsTerm(:, :, :, rxns(rxnIndex, 3)) = rhsTerm(:, :, :, rxns(rxnIndex, 3)) - individualReactionTerm;
        end

        % Second product (with first reactant only)
        if rxns(rxnIndex,4) ~= 0
            rhsTerm(:, :, :, rxns(rxnIndex, 4)) = rhsTerm(:, :, :, rxns(rxnIndex, 4)) - individualReactionTerm;
        end

    % Only second reactant is present
    elseif ((rxns(rxnIndex,1) == 0) && (rxns(rxnIndex,2) ~= 0))
        individualReactionTerm = -poro(1)*rxns(rxnIndex,5) ...
            * fields3DStar(nVar+rxns(rxnIndex,2):nVar:end-nVar, :, :);
        
        rhsTerm(:, :, :, rxns(rxnIndex, 2)) = rhsTerm(:, :, :, rxns(rxnIndex, 2)) + individualReactionTerm;

        % First product (with second reactant only)
        if rxns(rxnIndex,3) ~= 0
            rhsTerm(:, :, :, rxns(rxnIndex, 3)) = rhsTerm(:, :, :, rxns(rxnIndex, 3)) - individualReactionTerm;
        end

        % Second product (with second reactant only)
        if rxns(rxnIndex,4) ~= 0
            rhsTerm(:, :, :, rxns(rxnIndex, 4)) = rhsTerm(:, :, :, rxns(rxnIndex, 4)) - individualReactionTerm;
        end

    % Add constant term to RHS for constant production rate
    else
        individualReactionTerm = poro(1)*rxns(rxnIndex,5);
        % First product (no reactants)
        if rxns(rxnIndex,3) ~= 0
            rhsTerm(:, :, :, rxns(rxnIndex, 3)) = rhsTerm(:, :, :, rxns(rxnIndex, 3)) + individualReactionTerm;
        end

        % Second product (no reactants)
        if rxns(rxnIndex,4) ~= 0
            rhsTerm(:, :, :, rxns(rxnIndex, 4)) = rhsTerm(:, :, :, rxns(rxnIndex, 4)) + individualReactionTerm;
        end
    end
end

end