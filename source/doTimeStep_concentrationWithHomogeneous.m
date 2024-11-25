function [dxStripConcentrationOnly, t, residual] = doTimeStep_concentrationWithHomogeneous(dxC, dxF, nSpecies, xStripMiddle, xStripOld, xStripOldOld, xStripOldOldOld, isFrontNoFluxCondition, isBackNoFluxCondition,  ...
    t, dt, dtPreviousTimeStep, dtPrevPrevTimeStep, constants, faradaicRxnInfo, isNeutralGas, ...
    implicitTransportRHS, homogeneousReactionRHS, faradaicReactionRHS, implicitTransportRHSOld, explicitTransportRHSOld, homogeneousReactionRHSOld, faradaicReactionRHSOld, explicitTransportRHSOldOld, explicitTransportRHSOldOldOld, stepCtr)
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
vale = constants.vale;
e = constants.e;
kb = constants.kb;
T = constants.T;
NA = constants.NA;
phiElectrode = constants.phiElectrode;
nVar = nSpecies+1;

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

% Unsteady term computation
if stepCtr == 1
    timeAlpha = 1/dt;
    timeBeta = -1/dt;
    timeGamma = 0;
    timeDelta = 0;
    lhsMult = 1/2;
elseif stepCtr == 2
    timeRatio = dtPreviousTimeStep/dt;
    timeAlpha = 3/(2*dt);
    timeBeta = -3/2*(timeRatio + 1)^2/(timeRatio*(timeRatio + 2)) / dt;
    timeGamma = 3/2*1/(timeRatio*(timeRatio + 2)) / dt;
    timeDelta = 0;
    lhsMult = 1;
else
    timeRatio1 = dtPreviousTimeStep/dt;
    timeRatio2 = (dtPreviousTimeStep + dtPrevPrevTimeStep)/dt;
    timeAlpha = (11/6) / dt;
    timeBeta = -(11*timeRatio1^2*timeRatio2^2 + timeRatio1^2 + timeRatio1*timeRatio2 + 7*timeRatio1 + timeRatio2^2 + 7*timeRatio2)/(6*timeRatio1^2*timeRatio2^2) / dt;
    timeGamma = (timeRatio2 + 7)/(6*timeRatio1^2*(timeRatio2 - timeRatio1)) / dt;
    timeDelta = -(timeRatio1 + 7)/(6*timeRatio2^2*(timeRatio2 - timeRatio1)) / dt;
    lhsMult = 1;
end

% Determine LHS matrix A values, transport terms
for sIndex = 1:nSpecies
    mainDiag = zeros(xN-2, 1);
    subDiag = zeros(xN-3, 1);
    superDiag = zeros(xN-3, 1);

    effectiveDiff = poro(1)*diff(1,sIndex)/tort(1);

    mainDiag = mainDiag + poro(2:end-1) * timeAlpha;

    % Diffusion Flux
    if isFrontNoFluxCondition && ~isNeutralGas(sIndex)
        mainDiagFront = 0;
    else
        mainDiagFront = 2 * effectiveDiff./dxF(1);
%         mainDiagFront = 0;
    end
    
    if isBackNoFluxCondition && ~isNeutralGas(sIndex)
        mainDiagBack = 0;
    else
        mainDiagBack = 2 * effectiveDiff./dxF(end);
%         mainDiagBack = 0;
    end

    mainDiag = mainDiag + lhsMult * 1./dxC(2:end-1) .* ([effectiveDiff./dxF(2:end-1); mainDiagBack] + [mainDiagFront; effectiveDiff./dxF(2:end-1)]);
    subDiag = subDiag - lhsMult * effectiveDiff./dxF(2:end-1)./dxC(3:end-1);
    superDiag = superDiag - lhsMult * effectiveDiff./dxF(2:end-1)./dxC(2:end-2);

    % Electromigration flux (Need to refigure this for non-neutral species,
    % if some boundary should have a Dirichlet condition for charged ions.
    if isFrontNoFluxCondition && ~isNeutralGas(sIndex)
        mainDiagFront = 0;
    else
        mainDiagFront = 0;
    end
    
    if isBackNoFluxCondition && ~isNeutralGas(sIndex)
        mainDiagBack = 0;
    else
        mainDiagBack = 0;
    end

    mainDiag = mainDiag - lhsMult * 1/2 ./ dxC(2:end-1) * (vale(1, sIndex)*e/(kb*T)) .* ([dPhidX(2:end-1) .* effectiveDiff ; mainDiagBack] - [mainDiagFront; dPhidX(2:end-1).* effectiveDiff]);
    subDiag = subDiag + lhsMult * 1/2 ./ dxC(3:end-1) .* effectiveDiff * (vale(1, sIndex)*e/(kb*T)) .* dPhidX(2:end-1);
    superDiag = superDiag - lhsMult * 1/2 ./ dxC(2:end-2) .* effectiveDiff * (vale(1, sIndex)*e/(kb*T)) .* dPhidX(2:end-1);

    % Placement into sparse matrix vectors
    vals(sIndex:nSpecies:nSpecies*(xN-2)) = mainDiag;
    vals(nSpecies*(xN-2)+sIndex:nSpecies:nSpecies*(xN-2 + xN-3)) = subDiag;
    vals(nSpecies*(xN-2 + xN-3)+sIndex:nSpecies:nSpecies*(xN-2 + 2*(xN-3))) = superDiag;
end

% Determine LHS Matrix values, homogeneous reaction terms
startInd = nSpecies*(xN-2 + 2*(xN-3));
elemPerCell = nSpecies*nSpecies;
endInd = nSpecies*(xN-2 + 2*(xN-3)) + nSpecies*nSpecies*(xN-2);
rxns = constants.rxns{1};
for rxnIndex = 1:size(rxns, 1)
    if ((rxns(rxnIndex,1) ~= 0) && (rxns(rxnIndex,2) ~= 0))
        vals(startInd+(rxns(rxnIndex,1)-1)*nSpecies+rxns(rxnIndex,1):elemPerCell:endInd) ...
            = vals(startInd+(rxns(rxnIndex,1)-1)*nSpecies+rxns(rxnIndex,1):elemPerCell:endInd) ...
            + lhsMult * poro(1)*rxns(rxnIndex,5) * xStripMiddle(nVar+rxns(rxnIndex,2):nVar:end-nVar);

        vals(startInd+(rxns(rxnIndex,1)-1)*nSpecies+rxns(rxnIndex,2):elemPerCell:endInd) ...
            = vals(startInd+(rxns(rxnIndex,1)-1)*nSpecies+rxns(rxnIndex,2):elemPerCell:endInd) ...
            + lhsMult * poro(1)*rxns(rxnIndex,5) * xStripMiddle(nVar+rxns(rxnIndex,1):nVar:end-nVar);

        vals(startInd+(rxns(rxnIndex,2)-1)*nSpecies+rxns(rxnIndex,1):elemPerCell:endInd) ...
            = vals(startInd+(rxns(rxnIndex,2)-1)*nSpecies+rxns(rxnIndex,1):elemPerCell:endInd) ...
            + lhsMult * poro(1)*rxns(rxnIndex,5) * xStripMiddle(nVar+rxns(rxnIndex,2):nVar:end-nVar);

        vals(startInd+(rxns(rxnIndex,2)-1)*nSpecies+rxns(rxnIndex,2):elemPerCell:endInd) ...
            = vals(startInd+(rxns(rxnIndex,2)-1)*nSpecies+rxns(rxnIndex,2):elemPerCell:endInd) ...
            + lhsMult * poro(1)*rxns(rxnIndex,5) * xStripMiddle(nVar+rxns(rxnIndex,1):nVar:end-nVar);
        
        if rxns(rxnIndex,3) ~= 0
            vals(startInd+(rxns(rxnIndex,3)-1)*nSpecies+rxns(rxnIndex,1):elemPerCell:endInd) ...
            = vals(startInd+(rxns(rxnIndex,3)-1)*nSpecies+rxns(rxnIndex,1):elemPerCell:endInd) ...
            - lhsMult * poro(1)*rxns(rxnIndex,5) * xStripMiddle(nVar+rxns(rxnIndex,2):nVar:end-nVar);

            vals(startInd+(rxns(rxnIndex,3)-1)*nSpecies+rxns(rxnIndex,2):elemPerCell:endInd) ...
            = vals(startInd+(rxns(rxnIndex,3)-1)*nSpecies+rxns(rxnIndex,2):elemPerCell:endInd) ...
            - lhsMult * poro(1)*rxns(rxnIndex,5) * xStripMiddle(nVar+rxns(rxnIndex,1):nVar:end-nVar);
        end

        if rxns(rxnIndex,4) ~= 0
            vals(startInd+(rxns(rxnIndex,4)-1)*nSpecies+rxns(rxnIndex,1):elemPerCell:endInd) ...
            = vals(startInd+(rxns(rxnIndex,4)-1)*nSpecies+rxns(rxnIndex,1):elemPerCell:endInd) ...
            - lhsMult * poro(1)*rxns(rxnIndex,5) * xStripMiddle(nVar+rxns(rxnIndex,2):nVar:end-nVar);

            vals(startInd+(rxns(rxnIndex,4)-1)*nSpecies+rxns(rxnIndex,2):elemPerCell:endInd) ...
            = vals(startInd+(rxns(rxnIndex,4)-1)*nSpecies+rxns(rxnIndex,2):elemPerCell:endInd) ...
            - lhsMult * poro(1)*rxns(rxnIndex,5) * xStripMiddle(nVar+rxns(rxnIndex,1):nVar:end-nVar);
        end

    elseif ((rxns(rxnIndex,1) ~= 0) && (rxns(rxnIndex,2) == 0))
        vals(startInd+(rxns(rxnIndex,1)-1)*nSpecies+rxns(rxnIndex,1):elemPerCell:endInd) ...
            = vals(startInd+(rxns(rxnIndex,1)-1)*nSpecies+rxns(rxnIndex,1):elemPerCell:endInd) ...
            + lhsMult * poro(1)*rxns(rxnIndex,5);

        if rxns(rxnIndex,3) ~= 0
            vals(startInd+(rxns(rxnIndex,3)-1)*nSpecies+rxns(rxnIndex,1):elemPerCell:endInd) ...
            = vals(startInd+(rxns(rxnIndex,3)-1)*nSpecies+rxns(rxnIndex,1):elemPerCell:endInd) ...
            - lhsMult * poro(1)*rxns(rxnIndex,5);
        end

        if rxns(rxnIndex,4) ~= 0
            vals(startInd+(rxns(rxnIndex,4)-1)*nSpecies+rxns(rxnIndex,1):elemPerCell:endInd) ...
            = vals(startInd+(rxns(rxnIndex,4)-1)*nSpecies+rxns(rxnIndex,1):elemPerCell:endInd) ...
            - lhsMult * poro(1)*rxns(rxnIndex,5);
        end

    elseif ((rxns(rxnIndex,1) == 0) && (rxns(rxnIndex,2) ~= 0))
        vals(startInd+(rxns(rxnIndex,2)-1)*nSpecies+rxns(rxnIndex,2):elemPerCell:endInd) ...
            = vals(startInd+(rxns(rxnIndex,2)-1)*nSpecies+rxns(rxnIndex,2):elemPerCell:endInd) ...
            + lhsMult * poro(1)*rxns(rxnIndex,5);

        if rxns(rxnIndex,3) ~= 0
            vals(startInd+(rxns(rxnIndex,3)-1)*nSpecies+rxns(rxnIndex,2):elemPerCell:endInd) ...
            = vals(startInd+(rxns(rxnIndex,3)-1)*nSpecies+rxns(rxnIndex,2):elemPerCell:endInd) ...
            - lhsMult * poro(1)*rxns(rxnIndex,5);
        end

        if rxns(rxnIndex,4) ~= 0
            vals(startInd+(rxns(rxnIndex,4)-1)*nSpecies+rxns(rxnIndex,2):elemPerCell:endInd) ...
            = vals(startInd+(rxns(rxnIndex,4)-1)*nSpecies+rxns(rxnIndex,2):elemPerCell:endInd) ...
            - lhsMult * poro(1)*rxns(rxnIndex,5);
        end

    end
end
    
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
    b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) - poro(1)*(timeAlpha * xStripMiddle((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)) ...
        + timeBeta * xStripOld((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)) ...
        + timeGamma * xStripOldOld((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)) ...
        + timeDelta * xStripOldOldOld((nSpecies+1)+sIndex:(nSpecies+1):end-(nSpecies+1)));

    if stepCtr == 1
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(1/2*implicitTransportRHS(:, 1, 1, sIndex) + 1/2*implicitTransportRHSOld(:, 1, 1, sIndex), [xN-2,1]);
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(1/2*homogeneousReactionRHS(:, 1, 1, sIndex) + 1/2*homogeneousReactionRHSOld(:, 1, 1, sIndex), [xN-2,1]);
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(1/2*faradaicReactionRHS(:, 1, 1, sIndex) + 1/2*faradaicReactionRHSOld(:, 1, 1, sIndex), [xN-2,1]);
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(explicitTransportRHSOld(:, 1, 1, sIndex), [xN-2,1]);
    elseif stepCtr == 2
        rhsPrefactor = 3/2*(timeRatio+1)/(timeRatio+2);
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(rhsPrefactor * implicitTransportRHS(:, 1, 1, sIndex), [xN-2,1]);
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(rhsPrefactor * homogeneousReactionRHS(:, 1, 1, sIndex), [xN-2,1]);
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(rhsPrefactor * faradaicReactionRHS(:, 1, 1, sIndex), [xN-2,1]);
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(3/2*(timeRatio+1)^2/(timeRatio*(timeRatio + 2)) * explicitTransportRHSOld(:, 1, 1, sIndex) - 3/2*(timeRatio+1)/(timeRatio*(timeRatio + 2)) * explicitTransportRHSOldOld(:, 1, 1, sIndex), [xN-2,1]);
    else
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(implicitTransportRHS(:, 1, 1, sIndex), [xN-2,1]);
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(homogeneousReactionRHS(:, 1, 1, sIndex), [xN-2,1]);
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(faradaicReactionRHS(:, 1, 1, sIndex), [xN-2,1]);
        b(sIndex:nSpecies:end) = b(sIndex:nSpecies:end) + reshape(((timeRatio1+1)*(timeRatio2+1)/(timeRatio1*timeRatio2)) * explicitTransportRHSOld(:, 1, 1, sIndex) ...
            + (-(1+timeRatio2)/(timeRatio1*(timeRatio2-timeRatio1))) * explicitTransportRHSOldOld(:, 1, 1, sIndex) ...
            + ((1+timeRatio1)/(timeRatio2*(timeRatio2-timeRatio1))) * explicitTransportRHSOldOldOld(:, 1, 1, sIndex), [xN-2,1]);
    end
end


% chargeNorm = sqrt(sum(chargeField.^2))

%% Solve for the delta and update
% spparms('spumoni', 0)
warning('off', 'MATLAB:nearlySingularMatrix')
warning('off', 'MATLAB:singularMatrix')
if isreal(vals) && isreal(b) && ~any(isnan(vals)) && ~any(isnan(b))
    A = sparse(rowInd, colInd, vals, neq, neq, nnz);
%     L = ichol(A);
%     alpha = max(sum(abs(A),2)./diag(A))-2;
%     alpha = 1e3;
%     L1 = ichol(A, struct('type','ict','droptol',1e-4,'diagcomp',alpha));

%     options = struct("type","ilutp","droptol",1e-6);    
%     [L,U] = ilu(A, options);

%    [P,R,C] = equilibrate(A);
%    B = R*P*A*C;
%    d = R*P*b;

%    [y, ~] = bicgstab(B, d, 1e-6, 200);
%    dxStripConcentrationOnly = C*y;
     dxStripConcentrationOnly = A\b;
%     dxStripConcentrationOnly = bicgstab(A,b,[],1000);
%     L = ilu(A);
%     dxStripConcentrationOnly = bicg(A,b,1e-7,1000,L,L');
else
    dxStripConcentrationOnly = NaN((xN-2)*nSpecies, 1);
    residual = NaN((xN-2)*nSpecies, 1);
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
end
