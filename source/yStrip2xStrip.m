function fields3DRotated = yStrip2xStrip(fields3DIn, nX, nY, nZ, nSpecies)
    fields3DRotated = zeros(nX * (nSpecies+1), nY, nZ);
    for xIndex = 1:nX
        for yIndex = 1:nY
            for zIndex = 1:nZ
               for sIndex =1:nSpecies+1
                   fields3DRotated((xIndex-1)*(nSpecies+1) + sIndex, yIndex, zIndex) ...
                       = fields3DIn(xIndex, (yIndex-1)*(nSpecies+1) + sIndex, zIndex);
               end
            end 
        end 
    end
end