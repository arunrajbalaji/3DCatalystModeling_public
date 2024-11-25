function XYSlice = produce2DXYSlice(fields3D, zIndex, speciesIndex, nSpecies)
    % for speciesIndex = nSpecies+1, will extract electric potential
    XYSlice = squeeze(fields3D(speciesIndex:(nSpecies+1):end , :, zIndex));
end
