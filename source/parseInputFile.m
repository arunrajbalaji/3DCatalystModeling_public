function [layerInfo, uniqueSpeciesNames, latexNames, rxnInfo, faradaicRxnInfo, gasLayerInfo] = parseInputFile(filePath)
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
% Exceptions
ME_fileRead = MException('parseInputFile:fileNotRead',...
    'Input file not found, or not in JSON format.');
ME_fileParse = MException('parseInputFile:fileNotParsed',...
    'Input file not parsed successfully, check for errors.');

% try
    fileContents = fileread(filePath);
    rawOutput = jsondecode(fileContents);
% catch
%     throw(ME_fileRead)
% end
% 
% try
    layerInfo = struct([]);

    % Begin by packaging up layer physical/transport properties
    nLayers = length(rawOutput.layers);
    for ii = 1:nLayers
        layerInfo(ii).LID = rawOutput.layers(ii).layerID;
        layerInfo(ii).poro = rawOutput.layers(ii).porosity;
        layerInfo(ii).tort = rawOutput.layers(ii).tortuosity;
        layerInfo(ii).perm = rawOutput.layers(ii).permittivity;
        layerInfo(ii).xInterfaces = rawOutput.xInterfaceLocs;
        layerInfo(ii).yInterfaces = rawOutput.yInterfaceLocs;
        layerInfo(ii).zInterfaces = rawOutput.zInterfaceLocs;
        layerInfo(ii).dyMin = rawOutput.dyMin;
        layerInfo(ii).dyMax = rawOutput.dyMax;
        layerInfo(ii).yLRC = rawOutput.yLRC;
        layerInfo(ii).dzMin = rawOutput.dzMin;
        layerInfo(ii).dzMax = rawOutput.dzMax;
        layerInfo(ii).zLRC = rawOutput.zLRC;
        layerInfo(ii).xFrontNoFluxCondition = rawOutput.xFrontNoFluxCondition;
        layerInfo(ii).xBackNoFluxCondition = rawOutput.xBackNoFluxCondition;
        layerInfo(ii).xFrontPotentialBCValue = rawOutput.xFrontPotentialBCValue;
        layerInfo(ii).xBackPotentialBCValue = rawOutput.xBackPotentialBCValue;
        layerInfo(ii).yLeftNoFluxCondition = rawOutput.yLeftNoFluxCondition;
        layerInfo(ii).yRightNoFluxCondition = rawOutput.yRightNoFluxCondition;
        layerInfo(ii).yRightElectrolyteOverride = rawOutput.yRightElectrolyteOverride;
        layerInfo(ii).yRightBulkDirichlet = rawOutput.yRightBulkDirichlet;
        layerInfo(ii).yLeftPotentialBCValue = rawOutput.yLeftPotentialBCValue;
        layerInfo(ii).yRightPotentialBCValue = rawOutput.yRightPotentialBCValue;
        layerInfo(ii).zBottomNoFluxCondition = rawOutput.zBottomNoFluxCondition;
        layerInfo(ii).zTopNoFluxCondition = rawOutput.zTopNoFluxCondition;
        layerInfo(ii).zBottomPotentialBCValue = rawOutput.zBottomPotentialBCValue;
        layerInfo(ii).zTopPotentialBCValue = rawOutput.zTopPotentialBCValue;
        layerInfo(ii).T = rawOutput.T;
        layerInfo(ii).e = rawOutput.e;
        layerInfo(ii).kb = rawOutput.kb;
        layerInfo(ii).NA = rawOutput.NA;
        layerInfo(ii).m3ToLiter = rawOutput.m3ToLiter;
        layerInfo(ii).stericA = rawOutput.stericA;
        layerInfo(ii).backCharge = rawOutput.layers(ii).backCharge;
        layerInfo(ii).dxMax = rawOutput.layers(ii).dxMax;
        layerInfo(ii).dxMin = rawOutput.layers(ii).dxMin;
        layerInfo(ii).gridSymmetry = rawOutput.layers(ii).gridSymmetry;
        
        layerInfo(ii).timeStepOverride = rawOutput.timeStepOverride;
        layerInfo(ii).dt = rawOutput.dt;
        layerInfo(ii).dtOut = rawOutput.dtOut;
        layerInfo(ii).dtPlot = rawOutput.dtPlot;
        layerInfo(ii).tEnd = rawOutput.tEnd;
        layerInfo(ii).tStart = rawOutput.tStart;
        layerInfo(ii).reactionRampEndTime = rawOutput.reactionRampEndTime;
        layerInfo(ii).chargeThreshold = rawOutput.chargeThreshold;
        
        layerInfo(ii).homogeneousEq = rawOutput.homogeneousEq;
        layerInfo(ii).heq_tEnd = rawOutput.heq_tEnd;
        layerInfo(ii).heq_dt = rawOutput.heq_dt;
        layerInfo(ii).phiElectrode = rawOutput.phiElectrode;
        
        layerInfo(ii).newStart = rawOutput.newStart;
        layerInfo(ii).restartFileTime = rawOutput.restartFileTime;
        layerInfo(ii).debugMode = rawOutput.debugMode;

        layerInfo(ii).waterReactionIndices = rawOutput.waterReactionIndices;
        layerInfo(ii).galvanostaticMode = rawOutput.galvanostaticMode;
        layerInfo(ii).externalCapacitance = rawOutput.externalCapacitance;
        layerInfo(ii).targetCurrent = rawOutput.targetCurrent;
    end

    allSpeciesNames = [];

    % Determine the number and order of unique species
    for ii = 1:nLayers
        nSpecies = length(rawOutput.layers(ii).species);
        for jj = 1:nSpecies
            newSpecies = string(upper(rawOutput.layers(ii).species(jj).name));
            allSpeciesNames = [allSpeciesNames; newSpecies];
        end
    end

    uniqueSpeciesNames = unique(allSpeciesNames);
    nSpeciesTotal = length(uniqueSpeciesNames);
    
    % Fill in species-specific constants for each layer
    for ii = 1:nLayers
        nSpecies = length(rawOutput.layers(ii).species);
        
        % Logical array for determining which species are in the layer
        nameCellArray = {rawOutput.layers(ii).species(:).name};
        matchedLayers = matches(uniqueSpeciesNames,...
            nameCellArray, 'IgnoreCase', true);
        layerInfo(ii).species = matchedLayers;
        
        % Arrays for transferring specific constants
        % Default value zero for diffusvity (species not present)
        % Default value one for activity (species not present)
        layerInfo(ii).acti = ones(length(matchedLayers),1);
        layerInfo(ii).diff = zeros(length(matchedLayers),1);
        layerInfo(ii).vale = zeros(length(matchedLayers),1);
        layerInfo(ii).init = zeros(length(matchedLayers),1);
        layerInfo(ii).isNeutralGas = zeros(length(matchedLayers),1);
        layerInfo(ii).latexNames = cell(length(matchedLayers),1);
        
        % BC Handling (Dirichlet or Flux-based)
        layerInfo(ii).xFrontBCValues = zeros(length(matchedLayers),1);
        layerInfo(ii).xBackBCValues = zeros(length(matchedLayers),1);
        layerInfo(ii).yLeftBCValues(matchedLayers) = zeros(length(matchedLayers),1);
        layerInfo(ii).yRightBCValues(matchedLayers) = zeros(length(matchedLayers),1);
        layerInfo(ii).zBottomBCValues(matchedLayers) = zeros(length(matchedLayers),1);
        layerInfo(ii).zTopBCValues(matchedLayers) = zeros(length(matchedLayers),1);
        
        for jj = 1:nSpecies
        
            matchedLayers = matches(uniqueSpeciesNames,...
                nameCellArray{jj}, 'IgnoreCase', true);
        
            layerInfo(ii).acti(matchedLayers) = [rawOutput.layers(ii).species(jj).activityCoeff]';
        
            layerInfo(ii).diff(matchedLayers) = [rawOutput.layers(ii).species(jj).diffusionCoeff]';
            
            layerInfo(ii).vale(matchedLayers) = [rawOutput.layers(ii).species(jj).valence]';
            
            layerInfo(ii).init(matchedLayers) = [rawOutput.layers(ii).species(jj).initVal]';
            
            layerInfo(ii).isNeutralGas(matchedLayers) = [rawOutput.layers(ii).species(jj).isNeutralGas]';

            layerInfo(ii).latexNames{matchedLayers} = {rawOutput.layers(ii).species(jj).latexName}';
            
            % BC Values, for Dirichlet or Flux
            layerInfo(ii).xFrontBCValues(matchedLayers) = rawOutput.m3ToLiter * [rawOutput.xFrontBCValues(jj)]';
            layerInfo(ii).xBackBCValues(matchedLayers) = rawOutput.m3ToLiter * [rawOutput.xBackBCValues(jj)]';
            layerInfo(ii).yLeftBCValues(matchedLayers) = rawOutput.m3ToLiter * [rawOutput.yLeftBCValues(jj)]';
            layerInfo(ii).yRightBCValues(matchedLayers) = rawOutput.m3ToLiter * [rawOutput.yRightBCValues(jj)]';
            layerInfo(ii).zBottomBCValues(matchedLayers) = rawOutput.m3ToLiter * [rawOutput.zBottomBCValues(jj)]';
            layerInfo(ii).zTopBCValues(matchedLayers) = rawOutput.m3ToLiter * [rawOutput.zTopBCValues(jj)]';
        end
    end

    % Here, we assume only a single layer (as is done in basically all of
    % the rest of the code)
    latexNames = layerInfo(1).latexNames;
    
    % Fill in reaction information for each layer-for each layer: for each
    % species in the layer, for each reaction, include the reactants to
    % produce the species (and corresponding rate coefficient), or the
    % co-reactant to deplete the species (and the rate coefficient).
    % Forward and backward reactions must be entered separately.
    rxnInfo = struct([]);
    for ii = 1:nLayers
        nReacts = length(rawOutput.layers(ii).reactions);
        rxnInfo(ii).speciesRxns = zeros(nReacts, 3);
        for jj = 1:nReacts
            reactantIndex = find(matches(uniqueSpeciesNames,...
                rawOutput.layers(ii).reactions(jj).reactants{1}, 'IgnoreCase', true) == 1);
            if ~isempty(reactantIndex)
                rxnInfo(ii).speciesRxns(jj,1) = reactantIndex;
            else
                rxnInfo(ii).speciesRxns(jj,1) = 0;
            end
            
            if length(rawOutput.layers(ii).reactions(jj).reactants) < 2
                rxnInfo(ii).speciesRxns(jj,2) = 0;
            else
                reactantIndex = find(matches(uniqueSpeciesNames,...
                rawOutput.layers(ii).reactions(jj).reactants{2}, 'IgnoreCase', true) == 1);
                if ~isempty(reactantIndex)
                    rxnInfo(ii).speciesRxns(jj,2) = reactantIndex;
                else
                    rxnInfo(ii).speciesRxns(jj,2) = 0;
                end
            end
            
            productIndex = find(matches(uniqueSpeciesNames,...
                rawOutput.layers(ii).reactions(jj).products{1}, 'IgnoreCase', true) == 1);
            if ~isempty(productIndex)
                rxnInfo(ii).speciesRxns(jj,3) = productIndex;
            else
                rxnInfo(ii).speciesRxns(jj,3) = 0;
            end
            
            if length(rawOutput.layers(ii).reactions(jj).products) < 2
                rxnInfo(ii).speciesRxns(jj,4) = 0;
            else
                productIndex = find(matches(uniqueSpeciesNames,...
                rawOutput.layers(ii).reactions(jj).products{2}, 'IgnoreCase', true) == 1);
                if ~isempty(productIndex)
                    rxnInfo(ii).speciesRxns(jj,4) = productIndex;
                else
                    rxnInfo(ii).speciesRxns(jj,4) = 0;
                end
            end
            
            if rxnInfo(ii).speciesRxns(jj,1) ~= 0 && rxnInfo(ii).speciesRxns(jj,2) ~= 0
                rxnInfo(ii).speciesRxns(jj,5) = rawOutput.layers(ii).reactions(jj).forwardRateCoeff / layerInfo(1).m3ToLiter;
            elseif rxnInfo(ii).speciesRxns(jj,1) ~= 0 || rxnInfo(ii).speciesRxns(jj,2) ~= 0
                rxnInfo(ii).speciesRxns(jj,5) = rawOutput.layers(ii).reactions(jj).forwardRateCoeff;
            else
                rxnInfo(ii).speciesRxns(jj,5) = rawOutput.layers(ii).reactions(jj).forwardRateCoeff * layerInfo(ii).m3ToLiter;
            end
            rxnInfo(ii).speciesRxns(jj,6) = rawOutput.layers(ii).reactions(jj).wienCoeff;
            
        end
    end
    
    faradaicRxnInfo = struct([]);
    nFaradaicReacts = length(rawOutput.layers(1).faradaicReactions);
    for rxnIndex = 1:nFaradaicReacts
        reactantIndex = find(matches(uniqueSpeciesNames,...
                rawOutput.layers(1).faradaicReactions(rxnIndex).reactants, 'IgnoreCase', true) == 1);
        if ~isempty(reactantIndex)
            faradaicRxnInfo(rxnIndex).reactants = reactantIndex;
        else
            faradaicRxnInfo(rxnIndex).reactants = 0;
        end

        nProducts = length(rawOutput.layers(1).faradaicReactions(rxnIndex).products);
        faradaicRxnInfo(rxnIndex).products = zeros(nProducts,1);
        for tempProdIndex = 1:nProducts
            productIndex = find(matches(uniqueSpeciesNames,...
                rawOutput.layers(1).faradaicReactions(rxnIndex).products(tempProdIndex), 'IgnoreCase', true) == 1);
            if ~isempty(productIndex)
                faradaicRxnInfo(rxnIndex).products(tempProdIndex) = productIndex;
            else
                faradaicRxnInfo(rxnIndex).products(tempProdIndex) = 0;
            end
        end
        faradaicRxnInfo(rxnIndex).moleReactantsPerMoleElectron = rawOutput.layers(1).faradaicReactions(rxnIndex).moleReactantsPerMoleElectron;
        faradaicRxnInfo(rxnIndex).moleProductsPerMoleElectron = rawOutput.layers(1).faradaicReactions(rxnIndex).moleProductsPerMoleElectron;
        faradaicRxnInfo(rxnIndex).exchangeCurrentDensity = rawOutput.layers(1).faradaicReactions(rxnIndex).exchangeCurrentDensity;
        faradaicRxnInfo(rxnIndex).surfaceRoughness = rawOutput.layers(1).faradaicReactions(rxnIndex).surfaceRoughness;
        faradaicRxnInfo(rxnIndex).activationEnergy = rawOutput.layers(1).faradaicReactions(rxnIndex).activationEnergy;
        faradaicRxnInfo(rxnIndex).transferCoefficient = rawOutput.layers(1).faradaicReactions(rxnIndex).transferCoefficient;
        faradaicRxnInfo(rxnIndex).standardElectrodePotential = rawOutput.layers(1).faradaicReactions(rxnIndex).standardElectrodePotential;
        faradaicRxnInfo(rxnIndex).cRef = rawOutput.layers(1).faradaicReactions(rxnIndex).cRef* layerInfo(ii).m3ToLiter;
    end

    % Gas phase contants
    gasLayerInfo.gasFlowInlet = rawOutput.gasFlowInlet;
    % Note: this is corrected in main, to exclude the ghost cells
    gasLayerInfo.gasXYArea = rawOutput.spanwiseLength ...
        * ((rawOutput.gasDiffusionHeight + (1-rawOutput.gasFingerRatio)*(rawOutput.yInterfaceLocs(end) - rawOutput.yInterfaceLocs(1)))*rawOutput.gasPorosity ...
        + rawOutput.gasChannelRatio*rawOutput.gasChannelHeight);
    gasLayerInfo.gasXYAdvectionArea = rawOutput.spanwiseLength * (rawOutput.gasChannelRatio*rawOutput.gasChannelHeight);
    % Note: this is corrected in main, to exclude the ghost cells
    gasLayerInfo.nFingers = rawOutput.spanwiseLength * rawOutput.gasFingerRatio / (rawOutput.xInterfaceLocs(end) - rawOutput.xInterfaceLocs(1));
    gasLayerInfo.gasTortuosity = rawOutput.gasTortuosity;
    gasLayerInfo.spanwiseLength = rawOutput.spanwiseLength;
    gasLayerInfo.gasPressure = rawOutput.gasPressure;
    gasLayerInfo.gasInitialMoleFraction = rawOutput.gasInitialMoleFraction;
    gasLayerInfo.gasSpeciesNames = rawOutput.gasSpeciesNames;
    gasLayerInfo.gasSpeciesLatex = rawOutput.gasSpeciesLatex;
    gasLayerInfo.henrysLawConstant = rawOutput.henrysLawConstant ./ 1e2;
    gasLayerInfo.gasFingerRatio = rawOutput.gasFingerRatio;
    gasLayerInfo.gasDiffusionHeight = rawOutput.gasDiffusionHeight;
    gasLayerInfo.gasPorosity = rawOutput.gasPorosity;
    gasLayerInfo.gasChannelHeight = rawOutput.gasChannelHeight;
    gasLayerInfo.gasChannelRatio = rawOutput.gasChannelRatio;
    
    % Flags for convenient processing CASE 1: Is tracked in electrolyte
    gasLayerInfo.isTrackedInElectrolyte = zeros(size(gasLayerInfo.gasSpeciesNames));
    for sIndex = 1:length(gasLayerInfo.gasSpeciesNames)
        speciesIndexInElectrolyte = find(matches(uniqueSpeciesNames,...
                gasLayerInfo.gasSpeciesNames(sIndex), 'IgnoreCase', true) == 1);
        if ~isempty(speciesIndexInElectrolyte)
            gasLayerInfo.isTrackedInElectrolyte(sIndex) = speciesIndexInElectrolyte;
        end
    end
    
    % Flags for convenient processing CASE 2: Is PRODUCT of homogeneous reaction
    gasLayerInfo.isHomogeneousProduct = zeros(size(gasLayerInfo.gasSpeciesNames));
    for sIndex = 1:length(gasLayerInfo.gasSpeciesNames)
        for rIndex = 1:nReacts
            productIndex = find(matches(gasLayerInfo.gasSpeciesNames(sIndex),...
                rawOutput.layers(1).reactions(rIndex).products{1}, 'IgnoreCase', true) == 1);
            if ~isempty(productIndex)
                gasLayerInfo.isHomogeneousProduct(sIndex) = 1;
            end
            if length(rawOutput.layers(1).reactions(rIndex).products) > 1
                productIndex = find(matches(gasLayerInfo.gasSpeciesNames(sIndex),...
                    rawOutput.layers(1).reactions(rIndex).products{2}, 'IgnoreCase', true) == 1);
                if ~isempty(productIndex)
                    gasLayerInfo.isHomogeneousProduct(sIndex) = 1;
                end
            end
        end
    end
    
    % Flags for convenient processing CASE 3: Is PRODUCT of faradaic reaction
    gasLayerInfo.isFaradaicProduct = cell(size(gasLayerInfo.gasSpeciesNames));
    for sIndex = 1:length(gasLayerInfo.gasSpeciesNames)
        gasLayerInfo.isFaradaicProduct{sIndex} = zeros(nFaradaicReacts,1);
        for fIndex = 1:nFaradaicReacts
            productIndex = find(matches(gasLayerInfo.gasSpeciesNames(sIndex),...
                rawOutput.layers(1).faradaicReactions(fIndex).products(1), 'IgnoreCase', true) == 1);
            if ~isempty(productIndex)
                gasLayerInfo.isFaradaicProduct{sIndex}(fIndex) = 1;
            end
            if length(rawOutput.layers(1).faradaicReactions(fIndex).products) > 1
                productIndex = find(matches(gasLayerInfo.gasSpeciesNames(sIndex),...
                rawOutput.layers(1).faradaicReactions(fIndex).products(2), 'IgnoreCase', true) == 1);
                if ~isempty(productIndex)
                    gasLayerInfo.isFaradaicProduct{sIndex}(fIndex) = 2;
                end
            end
        end
    end
end
    
    
    % Package numerical information.
    
% catch
%   throw(ME_fileParse)
% end

%------------- END OF CODE --------------