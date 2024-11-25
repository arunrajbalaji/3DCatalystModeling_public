function doOutput(xCenter, yCenter, zCenter, uniqueSpecies, fields3D, gasFields, gasVelocity, gasSpeciesNames, t, folderName, bcPotentialOld)
% doOutput: Output simulation data at a certain time step, as a binary
% file.
%
% doOutput(xCenter, xFace, uniqueSpecies, z, t, folderName)
%
% Inputs:
%       uniqueSpecies   - ordered list of unique species names, from
%                         parsing function.
%       layerInfo       - layer info struct, as output by parsing function.
%       xCenter         - List of cell centers, where values are stored.
%                         Should not include an interface between layers.
%
% Outputs:
%       constants.poro            - 1D array of porosity
%       constants.tort            - 1D array of tortuosity
%       constants.perm            - 1D array of electric permittivity
%       constants.diff            - 2D array of species diffusivities
%       constants.acti            - 2D array of species tortuosities
%
% Example: 
%       Line 1 of example
%       Line 2 of example
%       Line 3 of example
%
% Other m-files required: none
% MAT-files required: none
%
% See also: parseInputFile.m
%
% Author: Arunraj Balaji
% Stanford University, Mani Group
% email: abalaji@stanford.edu
% Last revision: 25-June-2021
%------------- BEGIN CODE --------------

% Create file for writing
fileName = [folderName 'time_' num2str(t) '.bin'];
fileID = fopen(fileName, 'w');

% Create data structure for writing
% Structure outline:

% (z Value) (              xCenter             )
% (       ) (             first Species        )
% (       ) (            second Species        )
% (       ) (                ...               )
% (yCenter) (                ...               )
% (       ) (                ...               )
% (       ) (          last species            )
% (       ) (          elec. potential         )

% Slices in x-y plane are repeated, once for each z-value

nSpecies = length(uniqueSpecies);
dataOutputMatrix = zeros(length(zCenter)*(length(yCenter)*(nSpecies + 1) + 1), (length(xCenter) + 1));

for zIndex = 1:length(zCenter)
    dataOutputMatrix((zIndex-1)*(length(yCenter)*(nSpecies+1)+1) + 1, 1) = zCenter(zIndex);
    dataOutputMatrix((zIndex-1)*(length(yCenter)*(nSpecies+1)+1) + 1, 2:end) = xCenter';
    dataOutputMatrix((zIndex-1)*(length(yCenter)*(nSpecies+1)+1) + 2:(zIndex-1)*(length(yCenter)*(nSpecies+1)+1) + 1 + length(yCenter), 1) = yCenter;

    for ii = 1:nSpecies+1
        XYSlice = produce2DXYSlice(fields3D, zIndex, ii, nSpecies);
        dataOutputMatrix((zIndex-1)*(length(yCenter)*(nSpecies+1)+1) + 2 + (ii-1)*(length(yCenter)):(zIndex-1)*(length(yCenter)*(nSpecies+1)+1) + 1 + (ii)*(length(yCenter)), 2:end) = XYSlice';
    end
    
end

% Write data
fwrite(fileID, dataOutputMatrix, 'double');

% Close file
fclose(fileID);

% Create file for writing
fileName = [folderName 'gasFields_time_' num2str(t) '.bin'];
fileID = fopen(fileName, 'w');
fwrite(fileID, gasFields, 'double');
fclose(fileID);

fileName = [folderName 'gasVelocity_time_' num2str(t) '.bin'];
fileID = fopen(fileName, 'w');
fwrite(fileID, gasVelocity, 'double');
fclose(fileID);

fileName = [folderName 'bcPotential_time_' num2str(t) '.bin'];
fileID = fopen(fileName, 'w');
fwrite(fileID, bcPotentialOld, 'double');
fclose(fileID);

%------------- END OF CODE --------------