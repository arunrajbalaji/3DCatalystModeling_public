function doPlotGasVelocity(gasVelocity, cellCenters, figNum)
    plotStr = {'-o', '-x', '-d', '-s', '-.o', ':o', '-.s', '-.^'};
    
    figure(figNum)
    clf
    
    hold on
    plot(cellCenters/1e-2, gasVelocity/1e-2, '-o', 'Color', 'k')
    hold off
    
    set(gca, 'ticklabelinterpreter', 'latex', 'box', 'on')
    xlabel('$$z$$ (cm)', 'interpreter','latex')
    ylabel('$$\bar{u_z}$$ (cm/s)', 'interpreter', 'latex')
    box on
    set(gca, 'fontsize', 14)
    drawnow
end