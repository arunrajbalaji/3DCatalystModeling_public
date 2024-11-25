function doPlotChemEquilibrium(fields3D, cellCenters, xIndex, yIndex, zIndex, dimensionFlag, uniqueSpecies, reactionsToPlot, rxnInfo, m3ToLiter, figNum)
    cMap = linspecer(length(reactionsToPlot));
    plotStr = {'-o', '-x', '-d', '-s', '-.o', ':o', '-.s', ':d', '-.x'};
    
    nSpecies = length(uniqueSpecies);
    nReactionsToPlot = length(reactionsToPlot);
    
    reactionNames = cell(nReactionsToPlot,1);

    figure(figNum)
    clf

    for yIndex = 1:size(fields3D, 2)
        for zIndex = 1:size(fields3D, 3)
            for ii = 1:nSpecies
                fields3D(ii:(nSpecies+1):end, yIndex, zIndex) = fields3D(ii:(nSpecies+1):end, yIndex, zIndex) * m3ToLiter;
            end
        end
    end
    
    for reactionIndex = reactionsToPlot
        % Extract data first

        reactionNames{reactionIndex} = 'Reactants';
        
        eqConstant = ones(size(cellCenters));

        for structureIndex = 1:2
            speciesIndex = rxnInfo(2*reactionIndex - 1,structureIndex);
            if speciesIndex ~= 0
                
                if structureIndex == 1
                    reactionNames{reactionIndex} = [reactionNames{reactionIndex} ': ' uniqueSpecies{speciesIndex}];
                else
                    reactionNames{reactionIndex} = [reactionNames{reactionIndex} ' + ' uniqueSpecies{speciesIndex}];
                end

                if dimensionFlag == 1
                    XYSlice = produce2DXYSlice(fields3D, zIndex, speciesIndex, nSpecies);
                    dataToPlot = squeeze(XYSlice(:, yIndex));
                    xLabelText = '$$x$$ (m)';
                elseif dimensionFlag == 2
                    XYSlice = produce2DXYSlice(fields3D, zIndex, speciesIndex, nSpecies);
                    dataToPlot = squeeze(XYSlice(xIndex, :));
                    xLabelText = '$$y$$ (m)';
                elseif dimensionFlag == 3
                    YZSlice = produce2DYZSlice(fields3D, xIndex, speciesIndex, nSpecies);
                    dataToPlot = squeeze(YZSlice(yIndex, :));
                    xLabelText = '$$z$$ (m)';
                end
                eqConstant = eqConstant ./ dataToPlot;
            end
        end

        for structureIndex = 3:4
            speciesIndex = rxnInfo(2*reactionIndex - 1,structureIndex);
            if speciesIndex ~= 0
                if dimensionFlag == 1
                    XYSlice = produce2DXYSlice(fields3D, zIndex, speciesIndex, nSpecies);
                    dataToPlot = squeeze(XYSlice(:, yIndex));
                    xLabelText = '$$x$$ (m)';
                elseif dimensionFlag == 2
                    XYSlice = produce2DXYSlice(fields3D, zIndex, speciesIndex, nSpecies);
                    dataToPlot = squeeze(XYSlice(xIndex, :));
                    xLabelText = '$$y$$ (m)';
                elseif dimensionFlag == 3
                    YZSlice = produce2DYZSlice(fields3D, xIndex, speciesIndex, nSpecies);
                    dataToPlot = squeeze(YZSlice(yIndex, :));
                    xLabelText = '$$z$$ (m)';
                end
                eqConstant = eqConstant .* dataToPlot;
            end
        end
        
        equationEQConstant = rxnInfo(2*reactionIndex - 1,5)/rxnInfo(2*reactionIndex,5);

        % Then plot data
        hold on
        plot(cellCenters, eqConstant/equationEQConstant, plotStr{reactionIndex}, 'Color', cMap(reactionIndex, :))
        hold off
    end
    
    set(gca, 'ticklabelinterpreter', 'latex', 'box', 'on')
    xlabel(xLabelText, 'interpreter','latex')
    ylabel('Q/K$$_{\textrm{eq}}$$', 'interpreter', 'latex')
    legend(reactionNames, 'interpreter', 'latex', 'location', 'northeast')
    drawnow

end