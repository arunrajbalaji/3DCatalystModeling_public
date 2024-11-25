function YZSlice = produce2DYZSlice(fields3D, xIndex, speciesIndex, nSpecies)
    % for speciesIndex = nSpecies+1, will extract electric potential
    YZSlice = squeeze(fields3D((xIndex-1)*(nSpecies+1) + speciesIndex , :, :));
end
