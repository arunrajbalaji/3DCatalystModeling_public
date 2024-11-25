function doPlot3DBox(fields3D, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, uniqueSpeciesLatex, speciesIndex, plotType, figNum)
% outerFaces - will perform slices at faces of domain (x=0, y=0, z=max)
% Otherwise, will slice into center of x, y and z
nSpecies = length(uniqueSpecies);

figure(figNum)
clf

dataToPlot = permute(fields3D(nSpecies+1+speciesIndex:nSpecies+1:end-(nSpecies+1), :, :), [1, 3, 2]);
maxValue = max(dataToPlot(:));
exponent = log10(maxValue);
if exponent >= -1
    labelString = 'Concentration (M)';
    labelStringPotential = 'Potential $$\phi$$ (V)';
    dataToPlot = dataToPlot/1;
elseif exponent < -1 && exponent >= -4
    labelString = 'Concentration (mM)';
    labelStringPotential = 'Potential $$\phi$$ (mV)';
    dataToPlot = dataToPlot/1e-3;
elseif exponent < -4 && exponent >= -7
     labelString = 'Concentration ($$\mu$$M)';
     labelStringPotential = 'Potential $$\phi$$ ($$\mu$$V)';
     dataToPlot = dataToPlot/1e-6;
else
    labelString = 'Concentration (M)';
    labelStringPotential = 'Potential $$\phi$$ (V)';
    dataToPlot = dataToPlot/1;
end

xCenter = xCenter/1e-6;
yCenter = yCenter/1e-6;
zCenter = zCenter/1e-2;

xFace = xFace/1e-6;
yFace = yFace/1e-6;
zFace = zFace/1e-2;

[zz, xx, yy] = meshgrid(zCenter, xCenter(2:end-1), yCenter);

if plotType==1
    s = slice(zz, xx, yy, dataToPlot, [zCenter(end)], [xCenter(2)], [yCenter(1)]);
elseif plotType==0
    s = slice(zz, xx, yy, dataToPlot, [zCenter(end)/2], [xCenter(end)/2], [yCenter(end)/2]);
elseif plotType==2
    s = slice(zz, xx, yy, dataToPlot, [], [], [yCenter(2), yCenter(end)/2, yCenter(end-1)]);
elseif plotType==3
    s = slice(zz, xx, yy, dataToPlot, [zCenter(1), zCenter(end)/4, zCenter(end)/2, zCenter(end)*3/4, zCenter(end)], [], []);
elseif plotType==4
    s = slice(zz, xx, yy, dataToPlot, [], [], [yCenter(2), yCenter(end)/2, yCenter(end-1)]);
elseif plotType==5
    s = slice(zz, xx, yy, dataToPlot, [], [], [yCenter(2), yCenter(end)/2, yCenter(end-1)]);
end
set(s, 'edgecolor', 'none')
c = colorbar;
set(c, 'ticklabelinterpreter', 'latex')

if speciesIndex <= nSpecies
    c.Label.String= [uniqueSpeciesLatex{speciesIndex} ' ' labelString];
    c.Label.Interpreter = 'latex';
%     c.Label.FontSize = 24;
else
    c.Label.String = labelStringPotential;
    c.Label.Interpreter = 'latex';
%     c.Label.FontSize = 24;
end

xlabel('$$z$$ (cm)', 'interpreter', 'latex');
ylabel('$$x$$ ($$\mu$$m)', 'interpreter', 'latex');
zlabel('$$y$$ ($$\mu$$m)', 'interpreter', 'latex');

set(gca, 'ticklabelinterpreter', 'latex', 'box', 'on', 'fontsize', 14, 'TickDir', 'in')

pbaspect([0.4 0.1 0.2])
zlim([yFace(1) yFace(end)])
ylim([xFace(1) xFace(end)])

% scaleFactor = 1.05;
% outerPos = get(gca, 'OuterPosition');
% increasedWidth = (scaleFactor - 1) * outerPos(3);
% outerPos(3) = scaleFactor * outerPos(3);
% set(gca, 'OuterPosition', outerPos);
% 
% innerPos = get(gca, 'Position');
% innerPos(1) = innerPos(1) + increasedWidth/2;
% set(gca, 'Position', innerPos);

drawnow

end
