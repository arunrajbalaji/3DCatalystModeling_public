function doAllPlots(fields3D, gasFields, gasVelocity, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, uniqueSpeciesLatex, gasSpecies, gasSpeciesLatex, constants, gasLayerInfo)

nSpecies = length(uniqueSpecies);

% Change concentration units to mol/L, for plotting.
for yIndex = 1:size(fields3D, 2)
    for zIndex = 1:size(fields3D, 3)
        for ii = 1:nSpecies
            fields3D(ii:(nSpecies+1):end, yIndex, zIndex) = fields3D(ii:(nSpecies+1):end, yIndex, zIndex) / constants.m3ToLiter;
        end
    end
end

% Change units to partial pressure for gas fields
for sIndex = 1:length(gasSpecies)
    gasFields(sIndex:length(gasSpecies):end) = gasFields(sIndex:length(gasSpecies):end) ...
        * gasLayerInfo.gasPressure / 1e5;
end

%%% TURN BACK ON FOR STANDARD PLOTS
%doPlot3DBox(fields3D, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, 2, 3, 6)
%doPlot3DBox(fields3D, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, 1, 3, 27)
% doPlot3DBox(fields3D, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, uniqueSpeciesLatex, 1, 4, 1)
% doPlot3DBox(fields3D, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, uniqueSpeciesLatex, 2, 4, 2)
% doPlot3DBox(fields3D, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, uniqueSpeciesLatex, 3, 4, 3)
% doPlot3DBox(fields3D, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, uniqueSpeciesLatex, 4, 4, 4)
% doPlot3DBox(fields3D, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, uniqueSpeciesLatex, 5, 4, 5)
% doPlot3DBox(fields3D, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, 6, 4, 6)
% doPlot3DBox(fields3D, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, 7, 4, 7)
% doPlot3DBox(fields3D, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, 8, 3, 8)
doPlot1DProfiles(fields3D, xCenter, ceil(length(xCenter)/2), ceil(length(yCenter)/2), ceil(length(zCenter)/2), 1, uniqueSpecies, [1], 1)
doPlot1DProfiles(fields3D, xCenter, ceil(length(xCenter)/2), ceil(length(yCenter)/2), ceil(length(zCenter)/2), 1, uniqueSpecies, [2], 2)
doPlot1DProfiles(fields3D, xCenter, ceil(length(xCenter)/2), ceil(length(yCenter)/2), ceil(length(zCenter)/2), 1, uniqueSpecies, [3], 3)
doPlot1DProfiles(fields3D, xCenter, ceil(length(xCenter)/2), ceil(length(yCenter)/2), ceil(length(zCenter)/2), 1, uniqueSpecies, [4], 4)
doPlot1DProfiles(fields3D, xCenter, ceil(length(xCenter)/2), ceil(length(yCenter)/2), ceil(length(zCenter)/2), 1, uniqueSpecies, [5], 5)
% doPlot1DProfiles(fields3D, xCenter, ceil(length(xCenter)/2), ceil(length(yCenter)/2), ceil(length(zCenter)/2), 1, uniqueSpecies, [6], 6)
% doPlot1DProfiles(fields3D, xCenter, ceil(length(xCenter)/2), ceil(length(yCenter)/2), ceil(length(zCenter)/2), 1, uniqueSpecies, [7], 7)
% doPlot1DProfiles(fields3D, xCenter, ceil(length(xCenter)/2), ceil(length(yCenter)/2), ceil(length(zCenter)/2), 1, uniqueSpecies, [8], 8)

% doPlot1DProfiles(fields3D, yCenter, ceil(length(xCenter)/2), ceil(length(yCenter)/2), 1, 2, uniqueSpecies, [2], 2)
% doPlot1DProfiles(fields3D, yCenter, ceil(length(xCenter)/2), ceil(length(yCenter)/2), ceil(length(zCenter)/2), 2, uniqueSpecies, [1, 2, 3, 4, 5, 6, 7], 2)
% doPlot1DProfiles(fields3D, zCenter, ceil(length(xCenter)/2), ceil(length(yCenter)/2), ceil(length(zCenter)/2), 3, uniqueSpecies, [1, 2, 3, 4, 5, 6, 7], 3)
% doPlot1DProfiles(fields3D, yCenter, ceil(length(xCenter)/2), ceil(length(yCenter)/2), 1, 2, uniqueSpecies, [8], 4)
% doPlot1DProfiles(fields3D, zCenter, ceil(length(xCenter)/2), ceil(length(yCenter)/2), ceil(length(zCenter)/2), 3, uniqueSpecies, [8], 55)
% doPlot1DProfiles(fields3D, xCenter, ceil(length(xCenter)/2), ceil(length(yCenter)/2), ceil(length(zCenter)/2), 1, uniqueSpecies, [8], 58)
% doPlotChemEquilibrium(fields3D, xCenter, ceil(length(xCenter)/2), ceil(length(yCenter)/2), ceil(length(zCenter)/2), 1, uniqueSpecies, 1:8, constants.rxns{1}, constants.m3ToLiter, 9)

doPlotGasProfile(gasFields, gasLayerInfo.gasPressure, zCenter, gasSpecies, gasSpeciesLatex, [1 2 3], 11)
doPlotGasVelocity(gasVelocity, zFace, 15)
%   
% doPlotBulkCharge(fields3D, constants, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, 3, 99)
% doPlotBulkCharge(fields3D, constants, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, 4, 100)
  
% doPlot2DProfiles(fields3D, xCenter, yCenter, zCenter, ceil(length(xCenter)/2), ceil(length(yCenter)/2), ceil(length(zCenter)/2), 3, uniqueSpecies, [7], 12)


end