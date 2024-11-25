function XZSlice = produce2DXZSlice(fields3D, yIndex, speciesIndex, nSpecies)
    % for speciesIndex = nSpecies+1, will extract electric potential
    XZSlice = squeeze(fields3D(speciesIndex:(nSpecies+1):end , yIndex, :));
end
