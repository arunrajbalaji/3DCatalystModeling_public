function [deltaVals, t, residual] = doTimeStep_concentrationHomogeneousPre(eqConcValuesStar, eqConcValuesOld, eqConcValuesOldOld, eqConcValuesOldOldOld, ...
    t, dt, dtPreviousTimeStep, dtPrevPrevTimeStep, constants, homogeneousRxnSourceTermsStrip, homogeneousRxnSourceTermsStrip_init, stepCtr)
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
nSpecies = length(eqConcValuesStar);

t = t + dt;

% Initialization of variables
nnz = nSpecies * nSpecies;
neq = nSpecies;
b = zeros(neq,1);
vals = zeros(nnz, 1);

% Determining row index and column index vectors (sparse matrix)
rowInd = reshape(repmat((1:neq), [nSpecies, 1]), [nSpecies*nSpecies, 1]);
colInd = reshape(repmat((1:nSpecies)', [1, nSpecies]), [nSpecies*nSpecies, 1]);

% % Unsteady term computation
% if stepCtr == 1
%     lhsTime = poro(1)/dt;
%     lhsMult = 1/2;
% else%if stepCtr == 2
%     timeGamma = 1/(dtPreviousTimeStep * (1 + dtPreviousTimeStep/dt));
%     timeBeta = -timeGamma*(1 + dtPreviousTimeStep/dt)^2;
%     timeAlpha = -timeGamma - timeBeta;
%     lhsTime =  poro(1) * timeAlpha;
%     lhsMult = 1;
% % else
% %     delta_1 = dt;
% %     delta_2 = dt + dtPreviousTimeStep;
% %     delta_3 = dt + dtPreviousTimeStep + dtPrevPrevTimeStep;
% %     timeBeta = -delta_2*delta_3/(delta_1*(delta_1 - delta_2)*(delta_1 - delta_3));
% %     timeGamma = delta_3*delta_1/(delta_2*(delta_1 - delta_2)*(delta_2 - delta_3));
% %     timeDelta = delta_1*delta_2/(delta_3*(delta_1 - delta_3)*(delta_3 - delta_2));
% %     timeAlpha = -timeBeta - timeGamma - timeDelta;
% %     lhsTime =  poro(1) * timeAlpha;
% %     lhsMult = 1;
% end

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
lhsTime =  poro(1) * timeAlpha;

vals(1:(nSpecies+1):end) = vals(1:(nSpecies+1):end) + lhsTime;


% Determine LHS Matrix values, homogeneous reaction terms
rxns = constants.rxns{1};
for rxnIndex = 1:size(rxns, 1)
    if ((rxns(rxnIndex,1) ~= 0) && (rxns(rxnIndex,2) ~= 0))
        vals((rxns(rxnIndex,1)-1)*nSpecies+rxns(rxnIndex,1)) ...
            = vals((rxns(rxnIndex,1)-1)*nSpecies+rxns(rxnIndex,1)) ...
            + lhsMult * poro(1)*rxns(rxnIndex,5) * eqConcValuesStar(rxns(rxnIndex,2));

        vals((rxns(rxnIndex,1)-1)*nSpecies+rxns(rxnIndex,2)) ...
            = vals((rxns(rxnIndex,1)-1)*nSpecies+rxns(rxnIndex,2)) ...
            + lhsMult * poro(1)*rxns(rxnIndex,5) * eqConcValuesStar(rxns(rxnIndex,1));

        vals((rxns(rxnIndex,2)-1)*nSpecies+rxns(rxnIndex,1)) ...
            = vals((rxns(rxnIndex,2)-1)*nSpecies+rxns(rxnIndex,1)) ...
            + lhsMult * poro(1)*rxns(rxnIndex,5) * eqConcValuesStar(rxns(rxnIndex,2));

        vals((rxns(rxnIndex,2)-1)*nSpecies+rxns(rxnIndex,2)) ...
            = vals((rxns(rxnIndex,2)-1)*nSpecies+rxns(rxnIndex,2)) ...
            + lhsMult * poro(1)*rxns(rxnIndex,5) * eqConcValuesStar(rxns(rxnIndex,1));
        
        if rxns(rxnIndex,3) ~= 0
            vals((rxns(rxnIndex,3)-1)*nSpecies+rxns(rxnIndex,1)) ...
            = vals((rxns(rxnIndex,3)-1)*nSpecies+rxns(rxnIndex,1)) ...
            - lhsMult * poro(1)*rxns(rxnIndex,5) * eqConcValuesStar(rxns(rxnIndex,2));

            vals((rxns(rxnIndex,3)-1)*nSpecies+rxns(rxnIndex,2)) ...
            = vals((rxns(rxnIndex,3)-1)*nSpecies+rxns(rxnIndex,2)) ...
            - lhsMult * poro(1)*rxns(rxnIndex,5) * eqConcValuesStar(rxns(rxnIndex,1));
        end

        if rxns(rxnIndex,4) ~= 0
            vals((rxns(rxnIndex,4)-1)*nSpecies+rxns(rxnIndex,1)) ...
            = vals((rxns(rxnIndex,4)-1)*nSpecies+rxns(rxnIndex,1)) ...
            - lhsMult * poro(1)*rxns(rxnIndex,5) * eqConcValuesStar(rxns(rxnIndex,2));

            vals((rxns(rxnIndex,4)-1)*nSpecies+rxns(rxnIndex,2)) ...
            = vals((rxns(rxnIndex,4)-1)*nSpecies+rxns(rxnIndex,2)) ...
            - lhsMult * poro(1)*rxns(rxnIndex,5) * eqConcValuesStar(rxns(rxnIndex,1));
        end

    elseif ((rxns(rxnIndex,1) ~= 0) && (rxns(rxnIndex,2) == 0))
        vals((rxns(rxnIndex,1)-1)*nSpecies+rxns(rxnIndex,1)) ...
            = vals((rxns(rxnIndex,1)-1)*nSpecies+rxns(rxnIndex,1)) ...
            + lhsMult * poro(1)*rxns(rxnIndex,5);

        if rxns(rxnIndex,3) ~= 0
            vals((rxns(rxnIndex,3)-1)*nSpecies+rxns(rxnIndex,1)) ...
            = vals((rxns(rxnIndex,3)-1)*nSpecies+rxns(rxnIndex,1)) ...
            - lhsMult * poro(1)*rxns(rxnIndex,5);
        end

        if rxns(rxnIndex,4) ~= 0
            vals((rxns(rxnIndex,4)-1)*nSpecies+rxns(rxnIndex,1)) ...
            = vals((rxns(rxnIndex,4)-1)*nSpecies+rxns(rxnIndex,1)) ...
            - lhsMult * poro(1)*rxns(rxnIndex,5);
        end

    elseif ((rxns(rxnIndex,1) == 0) && (rxns(rxnIndex,2) ~= 0))
        vals((rxns(rxnIndex,2)-1)*nSpecies+rxns(rxnIndex,2)) ...
            = vals((rxns(rxnIndex,2)-1)*nSpecies+rxns(rxnIndex,2)) ...
            + lhsMult * poro(1)*rxns(rxnIndex,5);

        if rxns(rxnIndex,3) ~= 0
            vals((rxns(rxnIndex,3)-1)*nSpecies+rxns(rxnIndex,2)) ...
            = vals((rxns(rxnIndex,3)-1)*nSpecies+rxns(rxnIndex,2)) ...
            - lhsMult * poro(1)*rxns(rxnIndex,5);
        end

        if rxns(rxnIndex,4) ~= 0
            vals((rxns(rxnIndex,4)-1)*nSpecies+rxns(rxnIndex,2)) ...
            = vals((rxns(rxnIndex,4)-1)*nSpecies+rxns(rxnIndex,2)) ...
            - lhsMult * poro(1)*rxns(rxnIndex,5);
        end

    end
end

% Construct right side vector (computed in separate function and fed in)
b = b - poro(1)*(timeAlpha * eqConcValuesStar + timeBeta * eqConcValuesOld ...
    + timeGamma * eqConcValuesOldOld + timeDelta * eqConcValuesOldOldOld);

if stepCtr == 1
    b = b + 1/2*homogeneousRxnSourceTermsStrip + 1/2*homogeneousRxnSourceTermsStrip_init;
elseif stepCtr == 2
    rhsPrefactor = 3/2*(timeRatio+1)/(timeRatio+2);
    b = b + rhsPrefactor * homogeneousRxnSourceTermsStrip;
else
    b = b + homogeneousRxnSourceTermsStrip;
end

%% Solve for the delta and update
% spparms('spumoni', 0)
%     warning('off', 'MATLAB:nearlySingularMatrix')
if isreal(vals) && isreal(b) && ~any(isnan(vals)) && ~any(isnan(b))
    A = sparse(rowInd, colInd, vals, neq, neq, nnz);
    deltaVals = A\b;
%     dxStripConcentrationOnly = bicgstab(A,b,[],1000);
%     L = ilu(A);
%     dxStripConcentrationOnly = bicg(A,b,1e-7,1000,L,L');
else
    deltaVals = NaN(size(xStripMiddle));
    residual = NaN(size(xStripMiddle));
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
