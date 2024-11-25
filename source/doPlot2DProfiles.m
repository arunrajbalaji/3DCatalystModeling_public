function doPlot2DProfiles(fields3D, xCenter, yCenter, zCenter, xIndex, yIndex, zIndex, normalDirection, uniqueSpecies, speciesIndex, figNum)

nSpecies = length(uniqueSpecies);

figure(figNum)
clf

% Extract data first
if normalDirection == 3
    dataToPlot = produce2DXYSlice(fields3D(nSpecies+1+1:end-(nSpecies+1), :, :), zIndex, speciesIndex, nSpecies);
    
    [firstVar, secondVar] = meshgrid(xCenter(2:end-1), yCenter);
    
    xLabelText = '$$x$$ (m)';
    yLabelText = '$$y$$ (m)';
elseif normalDirection == 2
    dataToPlot = produce2DXZSlice(fields3D(nSpecies+1+1:end-(nSpecies+1), :, :), yIndex, speciesIndex, nSpecies);
    
    [firstVar, secondVar] = meshgrid(xCenter(2:end-1), zCenter);
    
    xLabelText = '$$x$$ (m)';
    yLabelText = '$$z$$ (m)';
elseif normalDirection == 1
    dataToPlot = produce2DYZSlice(fields3D, xIndex, speciesIndex, nSpecies);
    
    [firstVar, secondVar] = meshgrid(yCenter, zCenter);
    
    xLabelText = '$$y$$ (m)';
    yLabelText = '$$z$$ (m)';
end

% Then plot data

colormap pink

p = pcolor(firstVar, secondVar, dataToPlot');
set(p, 'edgecolor', 'none')

set(gca, 'ticklabelinterpreter', 'latex', 'box', 'on')
xlabel(xLabelText, 'interpreter','latex')
ylabel(yLabelText, 'interpreter','latex')

c = colorbar;
set(c, 'ticklabelinterpreter', 'latex')

if speciesIndex <= nSpecies
    c.Label.String= 'Concentration (M)';
    c.Label.Interpreter = 'latex';
else
    c.Label.String = 'Potential $$\Delta V$$ (V)';
    c.Label.Interpreter = 'latex';
end

drawnow
disp(num2str(max(max(dataToPlot))))

end