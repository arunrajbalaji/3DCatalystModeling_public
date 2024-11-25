function measuredCurrent = computeDeviceCurrentDensity(fields3D, constants, dxC, dyC, dyF, dzC, nSpecies, gasLayerInfo, hIndex, ohIndex, waterRxnIndices, potentialBCValue)
% Physical constants
poro = constants.poro;
tort = constants.tort;
diff = constants.diff;
vale = constants.vale;

e = constants.e;
kb = constants.kb;
T = constants.T;

NA = constants.NA;
nVar = nSpecies + 1;

nX = length(dxC);
nZ = length(dzC);

currentDensity = zeros(nX, nZ);

effectiveArea = gasLayerInfo.spanwiseLength ^ 2 * gasLayerInfo.gasFingerRatio;
deviceArea = gasLayerInfo.spanwiseLength ^ 2;
singleFingerArea = (sum(dxC(2:end-1)))*(sum(dzC));

% Species current density
for sIndex = 1:nSpecies
    effectiveDiff = diff(1, sIndex) * poro(1) / tort(1);
    
    % Extract data from fields3D
    xPlaneMiddle = fields3D(:, end, :);
    xPLaneLeft = fields3D(:, end-1, :);
    xPlaneLeftLeft = fields3D(:, end-2, :);
    xPlaneRight = zeros(size(fields3D(:, 1, :)));

     % Mesh spacing parameters
    dyFRight = dyF(end);
    dyFLeft = dyF(end);
    dyFLeftLeft = dyF(end-1);
    dyCCenter = dyC(end);

    if sIndex == hIndex
        K_eq = constants.rxns{1}(waterRxnIndices(1),5)/constants.rxns{1}(waterRxnIndices(2),5);
        ghostCellOH = -dyCCenter/2/dyFLeft*xPLaneLeft(ohIndex:nVar:end, 1, :) + (dyCCenter/2 + dyFLeft)/(dyFLeft) * xPlaneMiddle(ohIndex:nVar:end, 1, :);
        xPlaneRight(hIndex:nVar:end, 1, :) = K_eq ./ ghostCellOH;
    elseif vale(1, sIndex) >= 0
        xPlaneRight(sIndex:nVar:end) = constants.initVal(1,sIndex);
    end

   

    if vale(1, sIndex) < 0
        delta_1_2 = dyCCenter/2;
        delta_3_2 = dyCCenter/2 + dyFLeft;
        delta_5_2 = dyCCenter/2 + dyFLeft + dyFLeftLeft;
        beta = (-1/(delta_5_2 - delta_1_2)*(delta_1_2^2 - delta_5_2^2)) / (-(delta_5_2 - delta_3_2)/(delta_5_2 - delta_1_2) * (delta_1_2^2 - delta_5_2^2) + delta_3_2^2 - delta_5_2^2);
        alpha = 1/(delta_5_2 - delta_1_2) - beta * (delta_5_2 - delta_3_2)/(delta_5_2 - delta_1_2);
        gamma = -alpha - beta;

        yFluxRightDiff = - effectiveDiff * (alpha * xPlaneMiddle(sIndex:(nSpecies+1):end, 1, :) + beta * xPLaneLeft(sIndex:(nSpecies+1):end, 1, :) + gamma * xPlaneLeftLeft(sIndex:(nSpecies+1):end, 1, :));

        yFluxRightElectro = -vale(1,sIndex)*e/(kb * T)*effectiveDiff * (-dyCCenter/2/dyFLeft * xPLaneLeft(sIndex:(nSpecies+1):end, 1, :) + (dyCCenter/2 + dyFLeft)/(dyFLeft) * xPlaneMiddle(sIndex:(nSpecies+1):end, 1, :)) ...
            *2 .* (potentialBCValue - xPlaneMiddle(nSpecies+1 : nSpecies+1 : end, 1, :))/dyFRight;
    else
        yFluxRightDiff = - effectiveDiff * 2 * (xPlaneRight(sIndex:(nSpecies+1):end, 1, :) - xPlaneMiddle(sIndex:(nSpecies+1):end, 1, :)) ./ dyFRight;
        yFluxRightElectro = - effectiveDiff * vale(1,sIndex)*e/(kb*T) * xPlaneRight(sIndex:(nSpecies+1):end, 1, :) ...
            * 2 .* (potentialBCValue - xPlaneMiddle(nSpecies+1 : nSpecies+1 : end, 1, :)) ./ dyFRight;
    end

    totalRightFlux = yFluxRightDiff + yFluxRightElectro;

    currentDensity = currentDensity + reshape(totalRightFlux, [nX, nZ]) * vale(1, sIndex) * e * NA * (effectiveArea / singleFingerArea) / deviceArea;
end

measuredCurrent = sum(currentDensity(2:end-1, :) .* repmat(dxC(2:end-1), [1 nZ]) .* repmat(dzC', [nX-2 1]), [1 2]);
end