function [dxStripConcentrationOnly, t, residual, A] = doTimeStep_concentrationCouplingLinear(dxC, dxF, nSpecies, xStripMiddle, xStripOld, xStripOldOld, xStripOldOldOld, isFrontNoFluxCondition, isBackNoFluxCondition,  ...
    t, dt, dtPreviousTimeStep, dtPrevPrevTimeStep, constants, faradaicRxnInfo, isNeutralGas, faradaicRxnSourceTermsStrip, faradaicRxnSourceTermsStrip_init, transportRHSTermsStrip, transportRHSTermsStrip_init, stepCtr, debugMode)
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
nVar = nSpecies+1;

stericACubed = constants.stericACubed;
stericOnOffVec = constants.stericOnOffVec;

xN = length(dxC);

t = t + dt;

% Initialization of variables
nnz = nSpecies*(1+nSpecies)*(xN-2) + nSpecies*(2*(xN-3));
neq = nSpecies*(xN-2);
b = zeros(neq,1);
vals = zeros(nnz, 1);

% Determining row index and column index vectors (sparse matrix)
mainDiagRow = (1:1:neq)';
mainDiagCol = (1:1:neq)';

subDiagRow = (nSpecies+1:1:neq)';
subDiagCol = (1:1:neq-nSpecies)';

superDiagRow = (1:1:neq-nSpecies)';
superDiagCol = (nSpecies+1:1:neq)';

rxnBlockRow = reshape(repmat((1:neq), [nSpecies, 1]), [nSpecies*nSpecies*(xN-2), 1]);
rxnBlockCol = reshape(repmat(repmat((1:nSpecies)', [1, nSpecies]), [1, xN-2]) + reshape(repmat(nSpecies*(0:(xN-2-1)), [nSpecies, 1]), [1, (xN-2)*nSpecies]), [nSpecies*nSpecies*(xN-2), 1]);

rowInd = [mainDiagRow; subDiagRow; superDiagRow; rxnBlockRow];
colInd = [mainDiagCol; subDiagCol; superDiagCol; rxnBlockCol];

% Convenience variables
dPhidX = (xStripMiddle(2*nVar:nVar:end) - xStripMiddle(nVar:nVar:end-nVar))./dxF;

% Determine LHS matrix A values, transport terms
for sIndex = 1:nSpecies
    mainDiag = zeros(xN-2, 1);
    subDiag = zeros(xN-3, 1);
    superDiag = zeros(xN-3, 1);

    effectiveDiff = poro(1).*diff(1, sIndex)./tort(1);

    % Unsteady term computation
    if stepCtr == 1
        mainDiag = mainDiag + poro(2:end-1)/dt;
        lhsMult = 1/2;
    else%if stepCtr == 2
        timeGamma = 1/(dtPreviousTimeStep * (1 + dtPreviousTimeStep/dt));
        timeBeta = -timeGamma*(1 + dtPreviousTimeStep/dt)^2;
        timeAlpha = -timeGamma - timeBeta;
        mainDiag = mainDiag + poro(2:end-1) * timeAlpha;
        lhsMult = 1;
    %else
    %  delta_1 = dt;
    %  delta_2 = dt + dtPreviousTimeStep;
    %  delta_3 = dt + dtPreviousTimeStep + dtPrevPrevTimeStep;
    %  timeBeta = -delta_2*delta_3/(delta_1*(delta_1 - delta_2)*(delta_1 - delta_3));
    %  timeGamma = delta_3*delta_1/(delta_2*(delta_1 - delta_2)*(delta_2 - delta_3));
    %  timeDelta = delta_1*delta_2/(delta_3*(delta_1 - delta_3)*(delta_3 - delta_2));
    %  timeAlpha = -timeBeta - timeGamma - timeDelta;
    %  mainDiag = mainDiag + poro(2:end-1) * timeAlpha;
    %  lhsMult = 1;
    end

    % Diffusion Flux
    if isFrontNoFluxCondition && ~isNeutralGas(sIndex)
        mainDiagFront = 0;
    else
%         mainDiagFront = 2 * effectiveDiff./dxF(1);
        mainDiagFront = 0;
    end
    
    if isBackNoFluxCondition && ~isNeutralGas(sIndex)
        mainDiagBack = 0;
    else
%         mainDiagBack = 2 * effectiveDiff./dxF(end);
        mainDiagBack = 0;
    end

    mainDiag = mainDiag + lhsMult * 1./dxC(2:end-1) .* ([effectiveDiff./dxF(2:end-1); mainDiagBack] + [mainDiagFront; effectiveDiff./dxF(2:end-1)]);
    subDiag = subDiag - lhsMult * effectiveDiff./dxF(2:end-1)./dxC(3:end-1);
    superDiag = superDiag - lhsMult * effectiveDiff./dxF(2:end-1)./dxC(2:end-2);

    % Electromigration flux
    if isFrontNoFluxCondition && ~isNeutralGas(sIndex)
        mainDiagFront = 0;
    else
        mainDiagFront = 0;
    end
    
    if isBackNoFluxCondition && ~isNeutralGas(sIndex)
        mainDiagBack =      0;
    else
        mainDiagBack = 0;
    end

    mainDiag = mainDiag - lhsMult * 1/2 ./ dxC(2:end-1) * (vale(1, sIndex)*e/(kb*T)) .* ([dPhidX(2:end-1); mainDiagBack] - [mainDiagFront; dPhidX(2:end-1)]);
    subDiag = subDiag + lhsMult * 1/2 ./ dxC(2:end-2) * (vale(1, sIndex)*e/(kb*T)) .* dPhidX(2:end-1);
    superDiag = superDiag - lhsMult * 1/2 ./ dxC(2:end-2) * (vale(1, sIndex)*e/(kb*T)) .* dPhidX(2:end-1);

    % Placement into sparse matrix vectors
    vals(sIndex:nSpecies:nSpecies*(xN-2)) = mainDiag;
    vals(nSpecies*(xN-2)+sIndex:nSpecies:nSpecies*(xN-2 + xN-3)) = subDiag;
    vals(nSpecies*(xN-2 + xN-3)+sIndex:nSpecies:nSpecies*(xN-2 + 2*(xN-3))) = superDiag;
end
    
% Loop over Faradaic reactions and add to right and left side as
% necessary.
startInd = nSpecies*(xN-2 + 2*(xN-3));
elemPerCell = nSpecies*nSpecies;
endInd = nSpecies*(xN-2 + 2*(xN-3)) + nSpecies*nSpecies*(xN-2);
for fRxnIndex = 1:length(faradaicRxnInfo)
    if faradaicRxnInfo(fRxnIndex).reactants ~= 0
        mainDiag =  lhsMult * 1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleReactantsPerMoleElectron*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
            * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
            * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(2*nVar:nVar:end-nVar) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
            * 1/faradaicRxnInfo(fRxnIndex).cRef .* ones(xN-2, 1);
        vals(startInd+(faradaicRxnInfo(fRxnIndex).reactants-1)*nSpecies+faradaicRxnInfo(fRxnIndex).reactants:elemPerCell:endInd) = ...
            vals(startInd+(faradaicRxnInfo(fRxnIndex).reactants-1)*nSpecies+faradaicRxnInfo(fRxnIndex).reactants:elemPerCell:endInd) + mainDiag;
    end

    if faradaicRxnInfo(fRxnIndex).products(1) ~= 0
        if faradaicRxnInfo(fRxnIndex).reactants ~= 0
            mainDiag = -lhsMult * 1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(1)*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(2*nVar:nVar:end-nVar) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                * 1/faradaicRxnInfo(fRxnIndex).cRef .* ones(xN-2, 1);
            vals(startInd+(faradaicRxnInfo(fRxnIndex).products(1)-1)*nSpecies+faradaicRxnInfo(fRxnIndex).reactants:elemPerCell:endInd) = ...
            vals(startInd+(faradaicRxnInfo(fRxnIndex).products(1)-1)*nSpecies+faradaicRxnInfo(fRxnIndex).reactants:elemPerCell:endInd) + mainDiag;
        end
    end

    if faradaicRxnInfo(fRxnIndex).products(2) ~= 0
        if faradaicRxnInfo(fRxnIndex).reactants ~= 0
            mainDiag = -lhsMult * 1/(NA*e)*faradaicRxnInfo(fRxnIndex).moleProductsPerMoleElectron(2)*abs(faradaicRxnInfo(fRxnIndex).exchangeCurrentDensity) * faradaicRxnInfo(fRxnIndex).surfaceRoughness ...
                * exp(-faradaicRxnInfo(fRxnIndex).activationEnergy/(NA * kb * T)) ...
                * exp(-faradaicRxnInfo(fRxnIndex).transferCoefficient*e/(kb*T)*(phiElectrode - xStripMiddle(2*nVar:nVar:end-nVar) - faradaicRxnInfo(fRxnIndex).standardElectrodePotential)) ...
                * 1/faradaicRxnInfo(fRxnIndex).cRef .* ones(xN-2, 1);
            vals(startInd+(faradaicRxnInfo(fRxnIndex).products(2)-1)*nSpecies+faradaicRxnInfo(fRxnIndex).reactants:elemPerCell:endInd) = ...
            vals(startInd+(faradaicRxnInfo(fRxnIndex).products(2)-1)*nSpecies+faradaicRxnInfo(fRxnIndex).reactants:elemPerCell:endInd) + mainDiag;
        end
    end
end

for sIndex = 1:nSpecies
    if stepCtr == 1
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) - poro(1)/dt*(xStripMiddle((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)) - xStripOld((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)));
    else%if stepCtr == 2
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) - poro(1)*(timeAlpha*xStripMiddle((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)) ...
            + timeBeta*xStripOld((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)) ...
            + timeGamma*xStripOldOld((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)));
    %else
    %    b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) - poro(1)*(timeAlpha*xStripMiddle((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)) ...
    %        + timeBeta*xStripOld((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)) ...
    %        + timeGamma*xStripOldOld((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)) ...
    %        + timeDelta*xStripOldOldOld((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)));
    end

    if stepCtr == 1
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) ...
            + 1/2*reshape(faradaicRxnSourceTermsStrip(:, 1, 1, sIndex), [xN-2,1]) ...
            + 1/2*reshape(faradaicRxnSourceTermsStrip_init(:, 1, 1, sIndex), [xN-2,1]);

        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) ...
            + 1/2*reshape(transportRHSTermsStrip(:, 1, 1, sIndex), [xN-2,1]) ...
            + 1/2*reshape(transportRHSTermsStrip_init(:, 1, 1, sIndex), [xN-2,1]);
    else
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(faradaicRxnSourceTermsStrip(:, 1, 1, sIndex), [xN-2,1]);
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(transportRHSTermsStrip(:, 1, 1, sIndex), [xN-2,1]);
    end
end

%% Solve for the delta and update
% spparms('spumoni', 0)
%     warning('off', 'MATLAB:nearlySingularMatrix')
if isreal(vals) && isreal(b) && ~any(isnan(vals)) && ~any(isnan(b))
    A = sparse(rowInd, colInd, vals, neq, neq, nnz);
    dxStripConcentrationOnly = A\b;
%     dxStripConcentrationOnly = bicgstab(A,b,[],1000);
%     L = ilu(A);
%     dxStripConcentrationOnly = bicg(A,b,1e-7,1000,L,L');
else
    dxStripConcentrationOnly = NaN((xN-2)*nSpecies, 1);
    residual = NaN(size(xStripMiddle));
    A = nan(nSpecies*(xN-2), nSpecies*(xN-2));
    return
end

residual = -b;

% if debugMode
%         [~, maxIndex] = max(abs(dxStripConcentrationOnly));
%         maxDelta = dxStripConcentrationOnly(maxIndex);
% else
%         maxDelta = 0;
% end

%------------- END OF CODE --------------
