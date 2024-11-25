function newFolderPathList = generateDataFolders(filePath, taskIndex, jobNumber, parameterName, parameterValues)
% parseInputFile: Reads input file in JSON format, parses file, and
% generates two objects that store information. One stores general
% information about the simulation, and the other stores information about
% each layer specifically. The layer info is stored as layer objects inside
% an array.
%
% [layerInfo, uniqueSpeciesNames] = parseInputFile(filePath)
%
% Inputs:
%       filePath    - Full path to input file
%
% Outputs:
%       x_center    - array of ordered mesh cell centers
%       x_face      - array of ordered mesh cell faces      
%
% Example: 
%       Line 1 of example
%       Line 2 of example
%       Line 3 of example
%
% Other m-files required: none
% MAT-files required: none
%
% See also:
%
% Author: Arunraj Balaji
% Stanford University, Mani Group
% email: abalaji@stanford.edu
% Last revision: 16-June-2021
%------------- BEGIN CODE --------------

fileContents = fileread(filePath);
[folderPath,~,~] = fileparts(filePath);
rawOutput = jsondecode(fileContents);

numCases = length(parameterValues);

newFolderPathList = cell(numCases, 1);    

for cIndex = 1:numCases    
    newFolderName = [parameterName '_' num2str(parameterValues(cIndex))];
    newFolderPath = [folderPath '/output/' newFolderName];
    newFolderPathList{cIndex} = [newFolderPath '/'];
    
    if taskIndex == 1
        if (~exist(newFolderPath, 'dir'))
            mkdir(newFolderPath)
        end
        
        if strcmp(parameterName, 'yRightPotentialBCValue')
            eval(['rawOutput.' parameterName ' = parameterValues(cIndex);'])
            rawOutput.voltageValues = [-10 0:0.05:parameterValues(cIndex)];
            rawOutput.voltageSweepRate = (length(rawOutput.voltageValues) - 1)/rawOutput.parallelVoltageRampEndTime;
        elseif strcmp(parameterName, 'surfaceRoughness')
            for rIndex = 1:length(rawOutput.layers(1).faradaicReactions)
                eval(['rawOutput.layers(1).faradaicReactions(rIndex).' parameterName ' = parameterValues(cIndex);'])
            end
        elseif strcmp(parameterName, 'initCarbonate')
            rawOutput.layers(1).species(4).initVal = parameterValues(cIndex);
            rawOutput.layers(1).species(5).initVal = 0;
            rawOutput.layers(1).species(7).initVal = 2*parameterValues(cIndex);
            rawOutput.yRightBCValues(7) = 2*parameterValues(cIndex);
        elseif strcmp(parameterName, 'layers(1).porosity')
            rawOutput.layers(1).porosity = parameterValues(cIndex);
            rawOutput.layers(1).tortuosity = parameterValues(cIndex) ^ (-0.5);
        elseif strcmp(parameterName, 'xInterfaceLocs(2)')
            rawOutput.xInterfaceLocs(2) = parameterValues(cIndex);
            rawOutput.layers(1).dxMax = parameterValues(cIndex)/12.5;
        elseif strcmp(parameterName, 'zInterfaceLocs(2)')
            rawOutput.zInterfaceLocs(2) = parameterValues(cIndex);
            rawOutput.dzMax = parameterValues(cIndex)/40;
            rawOutput.dzMin = parameterValues(cIndex)/40;
        elseif strcmp(parameterName, 'yInterfaceLocs(2)')
            rawOutput.yInterfaceLocs(2) = parameterValues(cIndex);
            rawOutput.dyMax = parameterValues(cIndex)/15;
            rawOutput.dyMin = parameterValues(cIndex)/15;
        elseif strcmp(parameterName, 'fixedVolumeFingerHeight')
            rawOutput.yInterfaceLocs(2) = parameterValues(cIndex);
            rawOutput.dyMax = parameterValues(cIndex)/15;
            rawOutput.dyMin = parameterValues(cIndex)/15;
            rawOutput.gasFingerRatio = 0.4 * (75e-6 / parameterValues(cIndex));
        elseif strcmp(parameterName, 'fixedVolumeFingerThickness')
            rawOutput.xInterfaceLocs(2) = parameterValues(cIndex);
            rawOutput.layers(1).dxMax = parameterValues(cIndex)/12.5;
            rawOutput.gasFingerRatio = 0.4 * (5e-6 / parameterValues(cIndex));
        else
            eval(['rawOutput.' parameterName ' = parameterValues(cIndex);'])
        end
        
        outputText = jsonencode(rawOutput, 'PrettyPrint', true);
        
        fid = fopen([newFolderPath '/inputFile_' num2str(jobNumber) '.json'], 'w');
        fprintf(fid, outputText);
        fclose(fid);
        
        if rawOutput.newStart == 0
            dataList = dir([folderPath '/*time_' num2str(rawOutput.restartFileTime) '*.bin']);
            startFileNames = {dataList.name}';
            for dIndex = 1:length(startFileNames)
                copyfile([folderPath '/' startFileNames{dIndex}], [newFolderPath '/' startFileNames{dIndex}]);
            end
        end
    end
end
