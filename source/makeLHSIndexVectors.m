function [rowInd, colInd, nnz, nEq] = makeLHSIndexVectors(nSpecies, N, constants, faradaicRxnInfo)
% makeLHSIndexVectors: Make the LHS matrix (A) row and column index vectors
% in the correct order (cell locations close together, with potential at
% the end of each block).
%
%   [rowInd, colInd] = makeLHSIndexVectors(nSpecies)
%
% Inputs:
%
%       nSpecies        - Number of uique species
%       N               - Number of cell centers (unknown)
%
% Outputs:
%       rowInd          - LHS matrix A row indices, in correct order
%       colInd          - LHS matrix A column indices, in correct order
%
% Example: 
%       Line 1 of example
%       Line 2 of example
%       Line 3 of example
%
% Other m-files required: none
% MAT-files required: none
%
% See also: doTimeStep.m
%
% Author: Arunraj Balaji
% Stanford University, Mani Group
% email: abalaji@stanford.edu
% Last revision: 13-July-2021
%------------- BEGIN CODE --------------

% Data allocation
nnz = (N-2)*(nSpecies) + 2*nSpecies*(N-2-1) + nSpecies*(nSpecies)*(N-2) + 2*nSpecies*(nSpecies)*(N-2 - 1);
nEq = (N-2)*(nSpecies);     % Edge cells not computed, Dirichlet BC
rowInd = zeros(nnz,1);
colInd = zeros(nnz,1);

% Loop over iterations
%% Determine LHS matrix A
% Loop over species
ctr = 1;

% Species equations
for jj = 1:nSpecies
    % Loop over points, ignoring first and last cells (Dirichlet)
    % ColInd in matrix will be off by one compared to kk!!!
    for kk = 2:N-1
        % species cell to left
        if kk > 2
            rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + jj;
            colInd(ctr) = (kk-1 - 1)*(nSpecies)  + jj - (nSpecies);
            ctr = ctr + 1;
        end
        % species cell to right
        if kk < N-1
            rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + jj;
            colInd(ctr) = (kk-1 - 1)*(nSpecies)  + jj + (nSpecies);
            ctr = ctr + 1;
        end
        
        % species cell center
        rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + jj;
        colInd(ctr) = (kk-1 - 1)*(nSpecies)  + jj;
        ctr = ctr + 1;
        
        % Coupling to other species, for steric effects
        for ll=1:nSpecies
            % Cell to left
            if kk > 2
                rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + jj;
                colInd(ctr) = (kk-1 - 1 - 1)*(nSpecies)  + ll;
                ctr = ctr + 1;
            end
            
            % Cell to right
            if kk < N-1
                rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + jj;
                colInd(ctr) = (kk-1)*(nSpecies)  + ll;
                ctr = ctr + 1;
            end
            
            % Cell center
            rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + jj;
            colInd(ctr) = (kk-1 - 1)*(nSpecies)  + ll;
            ctr = ctr + 1;
        end
    end
end

for kk = 2:N-1
   % Add reaction terms for each cell
   rxns = constants.rxns{kk};
   
   for ll = 1:size(rxns, 1)
       % Categorized by configuration of reactants
       if ((rxns(ll,1) ~= 0) && (rxns(ll,2) ~= 0))
           rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,1);
           colInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,1);
           ctr = ctr + 1;

           rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,1);
           colInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,2);
           ctr = ctr + 1;

           rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,2);
           colInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,1);
           ctr = ctr + 1;

           rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,2);
           colInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,2);
           ctr = ctr + 1;

           if rxns(ll,3) ~= 0
               rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,3);
               colInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,1);
               ctr = ctr + 1;

               rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,3);
               colInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,2);
               ctr = ctr + 1;

           end

           if rxns(ll,4) ~= 0
               rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,4);
               colInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,1);
               ctr = ctr + 1;

               rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,4);
               colInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,2);
               ctr = ctr + 1;
           end
       
       elseif ((rxns(ll,1) ~= 0) && (rxns(ll,2) == 0))
           rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,1);
           colInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,1);
           ctr = ctr + 1;
           
           if rxns(ll,3) ~= 0
               rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,3);
               colInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,1);
               ctr = ctr + 1;
           end
           
           if rxns(ll,4) ~= 0
               rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,4);
               colInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,1);
               ctr = ctr + 1;
           end
           
       elseif ((rxns(ll,1) == 0) && (rxns(ll,2) ~= 0))
           rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,2);
           colInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,2);
           ctr = ctr + 1;
           
           if rxns(ll,3) ~= 0
               rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,3);
               colInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,2);
               ctr = ctr + 1;
           end
           
           if rxns(ll,4) ~= 0
               rowInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,4);
               colInd(ctr) = (kk-1 - 1)*(nSpecies)  + rxns(ll,2);
               ctr = ctr + 1;
           end
           
       else
           % Add constant term to RHS, no need to implement this here
           
       end
   end
   
   for fRxnIndex = 1:length(faradaicRxnInfo)
       if faradaicRxnInfo(fRxnIndex).reactants ~= 0
           rowInd(ctr) = (kk-2)*nSpecies + faradaicRxnInfo(fRxnIndex).reactants;
           colInd(ctr) = (kk-2)*nSpecies + faradaicRxnInfo(fRxnIndex).reactants;
           ctr = ctr + 1;
       end
       
       if faradaicRxnInfo(fRxnIndex).products(1) ~= 0
           if faradaicRxnInfo(fRxnIndex).reactants ~= 0
               rowInd(ctr) = (kk-2)*nSpecies + faradaicRxnInfo(fRxnIndex).products(1);
               colInd(ctr) = (kk-2)*nSpecies + faradaicRxnInfo(fRxnIndex).reactants;
               ctr = ctr + 1;
           end
       end
       
       if faradaicRxnInfo(fRxnIndex).products(2) ~= 0
           if faradaicRxnInfo(fRxnIndex).reactants ~= 0
               rowInd(ctr) = (kk-2)*nSpecies + faradaicRxnInfo(fRxnIndex).products(2);
               colInd(ctr) = (kk-2)*nSpecies + faradaicRxnInfo(fRxnIndex).reactants;
               ctr = ctr + 1;
           end
       end
   end
   
end

%------------- END OF CODE --------------