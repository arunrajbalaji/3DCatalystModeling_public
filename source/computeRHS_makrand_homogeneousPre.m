function homogeneousRxnSourceTerms = computeRHS_makrand_homogeneousPre(nSpecies, fields3DStar, constants)
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

homogeneousRxnSourceTerms = zeros(nSpecies, 1);

% Homogeneous reaction terms
rxns = constants.rxns{1};

for rxnIndex = 1:size(rxns, 1)
    if ((rxns(rxnIndex,1) ~= 0) && (rxns(rxnIndex,2) ~= 0))

        individualReactionTerm = -poro(1)*rxns(rxnIndex,5) ...
            * fields3DStar(rxns(rxnIndex,1)) ...
            .* fields3DStar(rxns(rxnIndex,2));

        homogeneousRxnSourceTerms(rxns(rxnIndex, 1)) = homogeneousRxnSourceTerms(rxns(rxnIndex, 1)) + individualReactionTerm;

        homogeneousRxnSourceTerms(rxns(rxnIndex, 2)) = homogeneousRxnSourceTerms(rxns(rxnIndex, 2)) + individualReactionTerm;
        
        if rxns(rxnIndex,3) ~= 0
            homogeneousRxnSourceTerms(rxns(rxnIndex, 3)) = homogeneousRxnSourceTerms(rxns(rxnIndex, 3)) - individualReactionTerm;
        end

        if rxns(rxnIndex,4) ~= 0
            homogeneousRxnSourceTerms(rxns(rxnIndex, 4)) = homogeneousRxnSourceTerms(rxns(rxnIndex, 4)) - individualReactionTerm;
        end

    elseif ((rxns(rxnIndex,1) ~= 0) && (rxns(rxnIndex,2) == 0))
        individualReactionTerm = -poro(1)*rxns(rxnIndex,5) ...
            * fields3DStar(rxns(rxnIndex,1));
        
        homogeneousRxnSourceTerms(rxns(rxnIndex, 1)) = homogeneousRxnSourceTerms(rxns(rxnIndex, 1)) + individualReactionTerm;

        if rxns(rxnIndex,3) ~= 0
            homogeneousRxnSourceTerms(rxns(rxnIndex, 3)) = homogeneousRxnSourceTerms(rxns(rxnIndex, 3)) - individualReactionTerm;
        end

        if rxns(rxnIndex,4) ~= 0
            homogeneousRxnSourceTerms(rxns(rxnIndex, 4)) = homogeneousRxnSourceTerms(rxns(rxnIndex, 4)) - individualReactionTerm;
        end

    elseif ((rxns(rxnIndex,1) == 0) && (rxns(rxnIndex,2) ~= 0))
        individualReactionTerm = -poro(1)*rxns(rxnIndex,5) ...
            * fields3DStar(rxns(rxnIndex,2));
        
        homogeneousRxnSourceTerms(rxns(rxnIndex, 2)) = homogeneousRxnSourceTerms(rxns(rxnIndex, 2)) + individualReactionTerm;

        if rxns(rxnIndex,3) ~= 0
            homogeneousRxnSourceTerms(rxns(rxnIndex, 3)) = homogeneousRxnSourceTerms(rxns(rxnIndex, 3)) - individualReactionTerm;
        end

        if rxns(rxnIndex,4) ~= 0
            homogeneousRxnSourceTerms(rxns(rxnIndex, 4)) = homogeneousRxnSourceTerms(rxns(rxnIndex, 4)) - individualReactionTerm;
        end
    else
        % Add constant term to RHS
        individualReactionTerm = poro(1)*rxns(rxnIndex,5);
        if rxns(rxnIndex,3) ~= 0
            homogeneousRxnSourceTerms(rxns(rxnIndex, 3)) = homogeneousRxnSourceTerms(rxns(rxnIndex, 3)) + individualReactionTerm;
        end

        if rxns(rxnIndex,4) ~= 0
            homogeneousRxnSourceTerms(rxns(rxnIndex, 4)) = homogeneousRxnSourceTerms(rxns(rxnIndex, 4)) + individualReactionTerm;
        end
    end

end
    
end

%------------- END OF CODE --------------
