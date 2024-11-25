function fields3DRotated = xStrip2yStrip(fields3DIn, nX, nY, nZ, nSpecies)
    fields3DRotated = zeros(nX, nY * (nSpecies+1), nZ);
    for xIndex = 1:nX
        for yIndex = 1:nY
            for zIndex = 1:nZ
               for sIndex =1:nSpecies+1
                   fields3DRotated(xIndex, (yIndex-1)*(nSpecies+1) + sIndex, zIndex) ...
                       = fields3DIn((xIndex-1)*(nSpecies+1) + sIndex, yIndex, zIndex);
               end
            end 
        end 
    end
end