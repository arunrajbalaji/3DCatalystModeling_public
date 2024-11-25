function [constants] = genPhysConstArrays(uniqueSpecies, layerInfo, rxnInfo, xCenter)
% parseInputFile: Reads input file in JSON format, parses file, and
% generates two objects that store information. One stores general
% information about the simulation, and the other stores information about
% each layer specifically. The layer info is stored as layer objects inside
% an array.
%
% [constants] = genPhysConstArrays(uniqueSpecies, layerInfo, xCenter)
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
% Unpack preliminary data
interfaceLocs = layerInfo(1).xInterfaces;
nSpecies = length(uniqueSpecies);
nLayers = length(interfaceLocs) - 1;

% Non-species-dependent properties (porosity, tortuosity, permeability
constants.poro = zeros(size(xCenter));
constants.tort = zeros(size(xCenter));
constants.perm = zeros(size(xCenter));
constants.backCharge = zeros(size(xCenter));

for ii = 1:nLayers
    constants.poro((xCenter > interfaceLocs(ii))&(xCenter < interfaceLocs(ii+1))) ...
        = layerInfo(ii).poro;
    constants.tort((xCenter > interfaceLocs(ii))&(xCenter < interfaceLocs(ii+1))) ...
        = layerInfo(ii).tort;
    constants.perm((xCenter > interfaceLocs(ii))&(xCenter < interfaceLocs(ii+1))) ...
        = layerInfo(ii).perm;
    constants.backCharge((xCenter > interfaceLocs(ii))&(xCenter < interfaceLocs(ii+1))) ...
        = layerInfo(ii).backCharge;
end

constants.e = layerInfo(1).e;
constants.T = layerInfo(1).T;
constants.kb = layerInfo(1).kb;
constants.NA = layerInfo(1).NA;
constants.m3ToLiter = layerInfo(1).m3ToLiter;
constants.stericA = layerInfo(1).stericA;
constants.stericACubed = constants.NA * constants.stericA^3;

waterPerm = 6.934e-10;
l_B = constants.e^2/(4*pi*waterPerm*constants.kb*constants.T);
constants.wienBeta = l_B*constants.e/(constants.kb*constants.T);

constants.phiElectrode = layerInfo(1).phiElectrode;

% Species-dependent properties (diffusivity and activity)
constants.diff = zeros(length(xCenter), nSpecies);
constants.acti = zeros(length(xCenter), nSpecies);
constants.vale = zeros(length(xCenter), nSpecies);
constants.initVal = zeros(length(xCenter), nSpecies);

for ii = 1:nLayers
    for jj = 1:nSpecies
        % If species is not present, then diffusivity should be zero
        constants.diff((xCenter > interfaceLocs(ii))&(xCenter < interfaceLocs(ii+1)),jj)...
            = layerInfo(ii).diff(jj);
        
        % If species is not present, then activity should be one (in LOG)
        constants.acti((xCenter > interfaceLocs(ii))&(xCenter < interfaceLocs(ii+1)),jj)...
            = layerInfo(ii).acti(jj);
        
        % If species is not present, then valence should be zero (but
        % doesn't matter, since diffusivity is already set to zero). Here,
        % will just set all values for all of space equal to the value for
        % the last time valence is specified. Valence should be specified
        % as same for each species in every layer.
        if layerInfo(ii).vale(jj) ~= 0
            constants.vale(:,jj) = layerInfo(ii).vale(jj);
        end
        
        constants.initVal((xCenter > interfaceLocs(ii))&(xCenter < interfaceLocs(ii+1)),jj)...
            = layerInfo(ii).init(jj) * constants.m3ToLiter;
    end
end

constants.stericOnOffVec = [constants.vale(1,:).*sign(constants.vale(1,:))./max(ones(size(constants.vale(1,:))), abs(constants.vale(1,:)))]';

% Reaction rate constants, formulated as a function of space (therefore, 2D
% array for storage)
constants.rxns = cell(length(xCenter),1);
for ii = 1:length(xCenter)
    layerNumber = max(find((xCenter(ii) > layerInfo(1).xInterfaces)));
    constants.rxns{ii} = rxnInfo(layerNumber).speciesRxns;
end

%------------- END OF CODE --------------