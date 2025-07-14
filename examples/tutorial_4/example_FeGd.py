import magnon

atoms, interactions = magnon.input.create_interacting_system('FeGd.poscar', 'FeGd_moments', 'FeGd_exchange', )
interactions = interactions.symmetrize(atoms)

special = {
    'G': [0,0,0],
    'H': [0.5,-0.5,0.5],
    'N': [0,0,0.5],
    'P': [0.25,0.25,0.25],
}

path = atoms.get_cell().bandpath(path='GHNGPH', npoints=240, special_points=special)

spectrum = magnon.MagnonSpectrum(atoms, interactions, ham_prefactor=1)

bstruct = spectrum.get_band_structure(path)
bstruct.plot(emin=0, emax=200, filename='FeGd_bands.png', ylabel='energies [meV]')

k_points = bstruct.path.cartesian_kpts()
bands = bstruct.energies

print(k_points[0], bands[0][0])

bstruct.write("FeGd_bands.json")

