function doPlotBulkCharge(fields3D, constants, xCenter, yCenter, zCenter, xFace, yFace, zFace, uniqueSpecies, plotType, figNum)
% outerFaces - will perform slices at faces of domain (x=0, y=0, z=max)
% Otherwise, will slice into center of x, y and z
nSpecies = length(uniqueSpecies);

figure(figNum)
clf

nX = length(xCenter);
nY = length(yCenter);
nZ = length(zCenter);

netCharge = zeros(nX-2, nY, nZ);

for sIndex = 1:nSpecies
    netCharge = netCharge + constants.vale(1, sIndex) * fields3D(nSpecies+1+sIndex:nSpecies+1:end-(nSpecies+1), :, :);
end

dataToPlot = permute(netCharge, [1, 3, 2]);
maxValue = max(dataToPlot(:));
exponent = log10(maxValue);
if exponent >= -1
    labelString = 'Net charge concentration (M)';
    dataToPlot = dataToPlot/1;
elseif exponent < -1 && exponent >= -4
    labelString = 'Net charge concentration (mM)';
    dataToPlot = dataToPlot/1e-3;
elseif exponent < -4 && exponent >= -7
     labelString = 'Net charge concentration ($$\mu$$M)';
     dataToPlot = dataToPlot/1e-6;
else
    labelString = 'Net charge concentration (M)';
    dataToPlot = dataToPlot/1;
end

xCenter = xCenter/1e-6;
yCenter = yCenter/1e-6;
zCenter = zCenter/1e-2;

xFace = xFace/1e-6;
yFace = yFace/1e-6;
zFace = zFace/1e-2;

dxC = (xFace(2:end) - xFace(1:end-1)) * 1e-6;
dyC = (yFace(2:end) - yFace(1:end-1)) * 1e-6;
dzC = (zFace(2:end) - zFace(1:end-1)) * 1e-2;

for xIndex =1:nX
    for yIndex = 1:nY
        for zIndex = 1:nZ
            dvC(xIndex, yIndex, zIndex) = dxC(xIndex) * dyC(yIndex) * dzC(zIndex);
        end
    end
end

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

c.Label.String= labelString;
c.Label.Interpreter = 'latex';

xlabel('$$z$$ (cm)', 'interpreter', 'latex');
ylabel('$$x$$ ($$\mu$$m)', 'interpreter', 'latex');
zlabel('$$y$$ ($$\mu$$m)', 'interpreter', 'latex');

set(gca, 'ticklabelinterpreter', 'latex', 'box', 'on', 'fontsize', 14, 'TickDir', 'in')
pbaspect([0.4 0.1 0.2])
zlim([yFace(1) yFace(end)])
ylim([xFace(1) xFace(end)])
drawnow

disp(['Maximum absolute charge: ' num2str(max(abs(netCharge(:))))])
disp(['Mean absolute charge: ' num2str(mean(abs(netCharge(:))))])
disp(['Total charge: ' num2str(constants.m3ToLiter * sum(netCharge .* dvC(2:end-1, :, :), [1 2 3]))])

end