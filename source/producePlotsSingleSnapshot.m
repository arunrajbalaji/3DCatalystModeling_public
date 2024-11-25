function producePlotsSingleSnapshot(folderName, time)
% Load convenient information from main_general output    
load([folderName 'constants.mat'], 'constants')
load([folderName 'uniqueSpecies.mat'], 'uniqueSpecies')
load([folderName 'meshSize.mat'], 'nX', 'nY', 'nZ', 'xFace', 'yFace', 'zFace')
load([folderName 'gasLayerInfo.mat'], 'gasLayerInfo')


% Read binary output file
[fields3D, xCenter, yCenter, zCenter] = ...
    readBinaryOutput([folderName 'time_' num2str(time) '.bin'], nX, nY, nZ, length(uniqueSpecies));

fileID = fopen([folderName 'gasFields_time_' num2str(time) '.bin']);
gasFields = fread(fileID, 'double');
fclose (fileID);

fileID = fopen([folderName 'gasVelocity_time_' num2str(time) '.bin']);
gasVelocity = fread(fileID, 'double');
fclose (fileID);

% Produce all plots for single snapshot
doAllPlots(fields3D, gasFields, gasVelocity, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, gasLayerInfo.gasSpeciesNames, constants)
end