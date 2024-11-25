function [xStripMiddle, t] = doTimeStep_concentrationOnly(dxC, dxF, dyC, dyFLeft, dyFRight, dzC, dzFBottom, dzFTop, ...
    nSpecies, xStripMiddle, xStripOld, xStripLeft, xStripRight, xStripBottom, xStripTop, ...
    isFrontFluxCondition, isBackFluxCondition, isLeftFluxCondition, isRightFluxCondition, isBottomFluxCondition, isTopFluxCondition, ...
    t, dt, constants, rowInd, colInd, nnz, nEq, xFrontBCValues, xBackBCValues, yLeftBCValues, yRightBCValues, zBottomBCValues, zTopBCValues, ...
    faradaicRxnInfo, isNeutralGas, yRightElectrolyteOverride, xStripLeftLeft, dyFLeftLeft, bulkElectrolyteL, bulkElectrolyteC, faradaicRxnSourceTermsStrip, transverseRHSTermsStrip)
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
acti = constants.acti;
vale = constants.vale;
e = constants.e;
kb = constants.kb;
T = constants.T;
NA = constants.NA;
phiElectrode = constants.phiElectrode;

rxns = constants.rxns{1};

stericACubed = constants.stericACubed;
stericOnOffVec = constants.stericOnOffVec;

xN = length(dxC);

effectiveDiffVec = repmat(poro(1)/tort(1)*diff(1,:)', [xN-2, 1]);
valenceVec = repmat(vale(1,:)', [xN-2, 1]);

vals = zeros(nnz,1);

t = t + dt;

% Must zero out b in advance (for reaction terms)
b = zeros(nEq,1);

% Add reaction terms as necessary
% Loop over grid points
for xLocIndex = 2:xN-1
    % Add reaction terms for each cell
    
    
    % Loop over Faradaic reactions and add to right and left side as
    % necessary.
    for fRxnIndex = 1:length(faradaicRxnInfo)
        if faradaicRxnInfo(fRxnIndex).reactants ~= 0
            vals(ctr) =  1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleReactantsPerMoleElectron*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                * 1/faradaicRxnInfo(fRxnIndex).cRef;
            ctr = ctr + 1;
        end
        
        if faradaicRxnInfo(fRxnIndex).products(1) ~= 0
            if faradaicRxnInfo(fRxnIndex).reactants ~= 0
                vals(ctr) = - 1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(1)*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                    * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                    * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                    * 1/faradaicRxnInfo(fRxnIndex).cRef;
                ctr = ctr + 1;
            end
        end
        
        if faradaicRxnInfo(fRxnIndex).products(2) ~= 0
            if faradaicRxnInfo(fRxnIndex).reactants ~= 0
                vals(ctr) = -1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(2)*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                    * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                    * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(xLocIndex*(nSpecies+1)) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                    * 1/faradaicRxnInfo(fRxnIndex).cRef;
                ctr = ctr + 1;
            end
        end
    end
end

for sIndex = 1:nSpecies
    b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(faradaicRxnSourceTermsStrip(:, 1, 1, sIndex), [xN-2,1]);
    b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(transverseRHSTermsStrip(:, 1, 1, sIndex), [xN-2,1]);
end

nMax = nSpecies*(xN-2);

% center diagonal, transport
unsteadyCenterTerm = repmat(poro(2:end-1)/dt, [nSpecies, 1]);
diffusionCenterTerm = effectiveDiffVec .* repmat(1./dxC(2:end-1).*(1./dxF(2:end) + 1./dxF(1:end-1)), [nSpecies, 1]);
if isFrontFlux && ~isNeutralGas(speciesIndex)
    diffusionCenterTerm(1:nSpecies, 1) = effectiveDiffVec(1:nSpecies, 1) .* repmat(1./dxC(2).*(1./dxF(2)), [nSpecies, 1]);
end
if isBackFlux && ~isNeutralGas(speciesIndex)
    diffusionCenterTerm(end-nSpecies:end, 1) = effectiveDiffVec(end-nSpecies:end, 1) .* repmat(1./dxC(end-1).*(1./dxF(end-1)), [nSpecies, 1]);
end
diffusionSuperDiag = -effectiveDiffVec(1:end-nSpecies,1) .* repmat(1./dxC(2:end-2)./dxF(2:end-1), [nSpecies, 1]);
diffusionSubDiag = -effectiveDiffVec(nSpecies+1:end,1) .* repmat(1./dxC(3:end-1)./dxF(2:end-1), [nSpecies, 1]);

electromigrationCenterTerm = -e/(kb*T) * valenceVec .* effectiveDiffVec .* repmat(1./dxC(2:end-1).*(1/2*dPhidX(2:end) - 1/2*dPhidX(1:end-1)), [nSpecies, 1]);
if isFrontFlux && ~isNeutralGas(speciesIndex)
    electromigrationCenterTerm(1:nSpecies, 1) = -e/(kb*T) * valenceVec(1:nSpecies,1) .* effectiveDiffVec(1:nSpecies,1) .* repmat(1./dxC(2).*(1/2*dPhidX(2)), [nSpecies, 1]);
end
if isBackFlux && ~isNeutralGas(speciesIndex)
    electromigrationCenterTerm(end-nSpecies:end, 1) = -e/(kb*T) * valenceVec(end-nSpecies:end,1) .* effectiveDiffVec(end-nSpecies:end,1) .* repmat(1./dxC(end-1).*(-1/2*dPhidX(end-1)), [nSpecies, 1]);
end
electromigrationSuperDiag = e/(kb*T) * valenceVec(1:end-nSpecies,1) .* effectiveDiffVec(1:end-nSpecies,1) .* repmat(1./dxC(2:end-2).*(1/2*dPhidX(2:end-1)), [nSpecies, 1]);
electromigrationSubDiag = e/(kb*T) * valenceVec(nSpecies+1:end,1) .* effectiveDiffVec(nSpecies+1:end,1) .* repmat(1./dxC(3:end-1).*(-1/2*dPhidX(2:end-1)), [nSpecies, 1]);

reactionDiagonals = zeros((xN-2)*nSpecies, nSpecies);

for rxnIndex = 1:size(rxns, 1)
    % Reactants side first
    if ((rxns(rxnIndex,1) ~= 0) && (rxns(rxnIndex,2) ~= 0))
        
        reactionDiagonals() = reactionDiagonals() + poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
        
        
        vals(ctr) = poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
        ctr = ctr + 1;
        
        vals(ctr) = poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
        ctr = ctr + 1;
        
        vals(ctr) = poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
        ctr = ctr + 1;
        
        vals(ctr) = poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
        ctr = ctr + 1;
        
        b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,1)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,1))...
            - poro(xLocIndex)*rxns(rxnIndex,5)...
            * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))...
            * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
        
        b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,2)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,2))...
            - poro(xLocIndex)*rxns(rxnIndex,5)...
            * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))...
            * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
        
        if rxns(rxnIndex,3) ~= 0
            vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
            ctr = ctr + 1;
            
            vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
            ctr = ctr + 1;
            
            b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3))...
                + poro(xLocIndex)*rxns(rxnIndex,5)...
                * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))...
                * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
        end
        
        if rxns(rxnIndex,4) ~= 0
            vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
            ctr = ctr + 1;
            
            vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5)*xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
            ctr = ctr + 1;
            
            b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4))...
                + poro(xLocIndex)*rxns(rxnIndex,5)...
                * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1))...
                * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
        end
        
    elseif ((rxns(rxnIndex,1) ~= 0) && (rxns(rxnIndex,2) == 0))
        vals(ctr) = poro(xLocIndex)*rxns(rxnIndex,5);
        ctr = ctr + 1;
        
        b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,1)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,1))...
            - poro(xLocIndex)*rxns(rxnIndex,5) * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
        
        
        if rxns(rxnIndex,3) ~= 0
            vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5);
            ctr = ctr + 1;
            
            b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3))...
                + poro(xLocIndex)*rxns(rxnIndex,5) * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
            
        end
        
        if rxns(rxnIndex,4) ~= 0
            vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5);
            ctr = ctr + 1;
            
            b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4))...
                + poro(xLocIndex)*rxns(rxnIndex,5) * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,1));
            
        end
        
    elseif ((rxns(rxnIndex,1) == 0) && (rxns(rxnIndex,2) ~= 0))
        vals(ctr) = poro(xLocIndex)*rxns(rxnIndex,5);
        ctr = ctr + 1;
        
        b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,2)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,2))...
            - poro(xLocIndex)*rxns(rxnIndex,5) * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
        
        if rxns(rxnIndex,3) ~= 0
            vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5);
            ctr = ctr + 1;
            
            b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,3))...
                + poro(xLocIndex)*rxns(rxnIndex,5) * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
        end
        
        if rxns(rxnIndex,4) ~= 0
            vals(ctr) = -poro(xLocIndex)*rxns(rxnIndex,5);
            ctr = ctr + 1;
            
            b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4)) = b((xLocIndex-1 - 1)*(nSpecies) + rxns(rxnIndex,4))...
                + poro(xLocIndex)*rxns(rxnIndex,5) * xStripMiddle((xLocIndex-1)*(nSpecies+1) + rxns(rxnIndex,2));
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


centerDiagonal = unsteadyCenterTerm + diffusionCenterTerm + electromigrationCenterTerm;
superDiagonal = diffusionSuperDiag + electromigrationSuperDiag;
subDiagonal = diffusionSubDiag + electromigrationSubDiag;


%% Solve for the delta and update
spparms('spumoni', 0)
%     warning('off', 'MATLAB:nearlySingularMatrix')

if isreal(vals) && isreal(b) && ~any(isnan(vals)) && ~any(isnan(b))
    A  = spdiags(diagonalEntries, [-fliplr(1:1:2*nSpecies-1) 0:1:2*nSpecies-1], nMax, nMax);
%     A = sparse(rowInd, colInd, vals, nEq, nEq, length(vals));
    dxStripConcentrationOnly = A\b;
else
    xStripMiddle = NaN(size(xStripMiddle));
    return
end

% Exit immediately if first iteration returns imaginary values
if ~isreal(dxStripConcentrationOnly) || any(isnan(dxStripConcentrationOnly))
    return
end


% Update the starred variable (concentration)
concentrationIndexingVector = [logical(zeros(nSpecies+1, 1)); ...
    repmat(logical([ones(nSpecies, 1); 0]), [xN-2, 1]); ...
    logical(zeros(nSpecies+1, 1))];
xStripMiddle(concentrationIndexingVector) ...
    = xStripMiddle(concentrationIndexingVector) + dxStripConcentrationOnly;

end

%------------- END OF CODE --------------
