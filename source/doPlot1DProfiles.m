function doPlot1DProfiles(fields3D, cellCenters, xIndex, yIndex, zIndex, dimensionFlag, uniqueSpecies, speciesToPlot, figNum)
    cMap = linspecer(length(speciesToPlot));
    plotStr = {'-o', '-x', '-d', '-s', '-.o', ':o', '-.s', '-.^'};
    
    nSpecies = length(uniqueSpecies);
    nSpeciesToPlot = length(speciesToPlot);
    
    figure(figNum)
    clf
    
    if dimensionFlag == 1
        cellCenters = cellCenters(2:end-1);
    end
    
    for speciesIndex = 1:nSpeciesToPlot
        % Extract data first
        if dimensionFlag == 1
            XYSlice = produce2DXYSlice(fields3D(nSpecies+1+1:end-(nSpecies+1), :, :), zIndex, speciesToPlot(speciesIndex), nSpecies);
            dataToPlot = squeeze(XYSlice(:, yIndex));
            xLabelText = '$$x$$ (m)';
        elseif dimensionFlag == 2
            XYSlice = produce2DXYSlice(fields3D, zIndex, speciesToPlot(speciesIndex), nSpecies);
            dataToPlot = squeeze(XYSlice(xIndex, :));
            xLabelText = '$$y$$ (m)';
        elseif dimensionFlag == 3
            YZSlice = produce2DYZSlice(fields3D, xIndex, speciesToPlot(speciesIndex), nSpecies);
            dataToPlot = squeeze(YZSlice(yIndex, :));
            xLabelText = '$$z$$ (m)';
        end
        
        % Then plot data
        hold on
        plot(cellCenters, dataToPlot, plotStr{speciesIndex}, 'Color', cMap(speciesIndex, :))
        hold off
    end
    
    set(gca, 'ticklabelinterpreter', 'latex', 'box', 'on')
    xlabel(xLabelText, 'interpreter','latex')
    
    if max(speciesToPlot) <= nSpecies
        ylabel('Concentration (M)', 'interpreter', 'latex')
        legend(uniqueSpecies(speciesToPlot), 'interpreter', 'latex', 'location', 'southwest')
    else
        ylabel('Potential $$\Delta V$$ (V)', 'interpreter', 'latex')
    end
    
    box on
    set(gca, 'fontsize', 14)
    drawnow

end