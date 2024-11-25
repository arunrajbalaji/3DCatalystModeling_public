function doPlotGasProfile(gasFields, gasPressure, cellCenters, gasSpecies, gasSpeciesLatex, speciesToPlot, figNum)
    cMap = linspecer(length(speciesToPlot));
    plotStr = {'-o', '-x', '-d', '-s', '-.o', ':o', '-.s', '-.^'};
    
    nSpecies = length(gasSpecies);
    nSpeciesToPlot = length(speciesToPlot);
    
    figure(figNum)
    clf
    
    for speciesIndex = 1:nSpeciesToPlot
        % Then plot data
        hold on
        plot(cellCenters/1e-2, gasFields(speciesIndex:nSpecies:end), ...
            plotStr{speciesIndex}, 'Color', cMap(speciesIndex, :), 'MarkerFaceColor', cMap(speciesIndex, :), 'linewidth', 1.5)
        hold off
    end
    
    set(gca, 'ticklabelinterpreter', 'latex', 'box', 'on')
    xlabel('$$z$$ (cm)', 'interpreter','latex')
    ylabel('$$p_i$$ (bar)', 'interpreter', 'latex')
    legend(gasSpeciesLatex{speciesToPlot}, 'interpreter', 'latex', 'location', 'northeast')
    box on
    set(gca, 'fontsize', 14)
    drawnow

end