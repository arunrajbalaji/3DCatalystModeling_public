function [fields3D, xCenter, yCenter, zCenter] = readBinaryOutput(fileName, nX, nY, nZ, nSpecies)
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

fileID = fopen(fileName);

rawData = fread(fileID, 'double');
fclose (fileID);

rawData = reshape(rawData, [nZ*(nY*(nSpecies + 1) + 1), (nX + 1)]);

xCenter = rawData(1, 2:end)';
yCenter = rawData(2:nY+1, 1);
zCenter = rawData(1:(nY*(nSpecies + 1) + 1):end, 1);

fields3D = zeros(nX*(nSpecies+1), nY, nZ);
for zIndex = 1:nZ
    for sIndex = 1:nSpecies+1
        XYSlice = transpose(rawData(((zIndex-1)*(nY*(nSpecies+1) + 1) + 1 + (sIndex-1)*nY + 1):((zIndex-1)*(nY*(nSpecies+1) + 1) + 1 + (sIndex-1)*nY + nY), 2:end));
        fields3D(sIndex:(nSpecies+1):end, :, zIndex) = XYSlice;
    end
end

end
%------------- END OF CODE --------------