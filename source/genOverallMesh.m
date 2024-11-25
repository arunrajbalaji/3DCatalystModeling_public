function [xCenter, xFace, dxC, dxF, yCenter, yFace, dyC, dyF, zCenter, zFace, dzC, dzF] = genOverallMesh(xLocs, dxMaxs, dxMins, LIDs, LRCs, yLocs, dyMin, dyMax, yLRC, zLocs, dzMin, dzMax, zLRC)
% genOverallMesh: Generates 1D non-uniform mesh points for a multi-layer
% geometry, by calling on appropriate individual layer mesh generation
% functions. Reads parameters that are provided as arrays, with an entry
% for each layer.
%
% Syntax:  [x_center, x_face] = genOverallMesh(xLocs, dxMaxs, dxMins, LIDs, LRCs)
%
% Inputs:
%       xLocs       - x - Locations of interfaces in increasing order,
%                     including domain boundaries. Size one larger than
%                     number of layers.
%       dxMaxs      - Maximum allowable mesh spacing (array, entry per layer)
%       dxMins      - Minimum allowable mesh spacing (array, entry per layer)
%       LRCs        - Left, right, center nonuniformity (array, entry per
%                     layer)
%
% Outputs:
%       x_center    - array of ordered mesh cell centers
%       x_face      - array of ordered mesh cell faces      
%       dx_c        - cell center grid sizes
%       dx_f        - cell face grid sizes
%
% Example: 
%       Line 1 of example
%       Line 2 of example
%       Line 3 of example
%
% Other m-files required: genLayerMesh.m
% MAT-files required: none
%
% See also: 
%
% Author: Arunraj Balaji
% Stanford University, Mani Group
% email: abalaji@stanford.edu
% Last revision: 16-June-2021
%------------- BEGIN CODE --------------

% Input checking and error throwing
ME_input = MException('genOverallMesh:incorrectInput', 'Input lengths are not aligned.');
if length(dxMaxs) ~= length(dxMins)
    throw(ME_input)
end

if length(dxMaxs) ~= length(LIDs)
    throw(ME_input)
end

if length(dxMaxs) ~= length(LRCs)
    throw(ME_input)
end

if isempty(dxMaxs)
    throw(ME_input)
end

if length(xLocs) ~= (length(dxMaxs)+1)
    throw(ME_input)
end

% Allocate vectors for output
xCenter = [];
xFace = [];

% Form input vectors for individual layer mesh generation function
xStarts = xLocs(1:end-1);
xEnds = xLocs(2:end);

%try
    
    for ii = 1:(length(dxMaxs))
        
        % Perform layer mesh generation
        [xCenterRaw, xFaceRaw] ...
            = genLayerMesh(xStarts(ii), xEnds(ii), dxMaxs(ii), dxMins(ii), LIDs(ii), LRCs(ii));
        
        % Assign centers
        xCenter = [xCenter; xCenterRaw];
        
        % Ignore last value each time, since same as first of next segment.
        xFace = [xFace; xFaceRaw(1:end-1)];
    end
    
    % Add last value
    xFace = [xFace; xEnds(end)];
    
    dxC = xFace(2:end) - xFace(1:end-1);
    dxF = xCenter(2:end) - xCenter(1:end-1);
    
    % Y-direction mesh generation starts here
    
    % Perform layer mesh generation
    [yCenter, yFace] ...
        = genLayerMesh(yLocs(1), yLocs(2), dyMax, dyMin, -1, yLRC);
    
    dyC = yFace(2:end) - yFace(1:end-1);
    dyF = yCenter(2:end) - yCenter(1:end-1);
    
    dyF = [dyF(1); dyF; dyF(end)];
    
    % Z-direction mesh generation starts here
    
    % Perform layer mesh generation
    [zCenter, zFace] ...
        = genLayerMesh(zLocs(1), zLocs(2), dzMax, dzMin, -1, zLRC);
    
    dzC = zFace(2:end) - zFace(1:end-1);
    dzF = zCenter(2:end) - zCenter(1:end-1);

    dzF = [dzF(1); dzF; dzF(end)];
    
%catch
%    warning('Error in genOverallMesh, assigning uniform mesh with N=100')
%    
%    x_face = (linspace(xStarts(1), xEnds(end), 100))';
%    x_center = 0.5*(x_face(2:end) + x_face(1:end-1));

%end

%------------- END OF CODE --------------
