import json
import logging
import unittest

import matplotlib.pyplot as plt
import numpy as np

import magnon

from tests.test_utils import band_structures_are_equal

logger = logging.getLogger("magnon")

class TestCrPS4(unittest.TestCase):
    """
    Uses exchange data from Phys. Rev. B 102, 024408 (2020) (https://dx.doi.org/10.1103/physrevb.102.024408) which
    also contains neutron scattering data to compare the magnon spectrum against.
    """
    do_plot = False
    do_overwrite = False

    logger.setLevel(logging.DEBUG)

    # Console handler for screen output
    console_handler = logging.StreamHandler()
    logger.addHandler(console_handler)

    def test_conventional_cell_full_exchange_spectrum(self):
        atoms, interactions = magnon.input.create_interacting_system(
            'tests/CrPS4/CrPS4.poscar', 'tests/CrPS4/CrPS4.moments',
            'tests/CrPS4/CrPS4.exchange', 1)
        interactions = interactions.symmetrize(atoms)

        special_kpts = {
            "G": [0, 0, 0],
            "X": [1, 0, 0],
            "V1": [1, 1, 0],
            "L1": [1, 1, 1],
            "L2": [0.5, 0.5, 0.5],
            "V2": [0.5, 0.5, 0.0]
        }

        path = atoms.get_cell().bandpath(special_kpts, npoints=50, special_points=special_kpts)

        magnon_spectrum = magnon.MagnonSpectrum(atoms, interactions, num_threads=2)
        bstruct = magnon_spectrum.get_band_structure(path)

        if self.do_overwrite:
            data = {
                "kpts": bstruct.path.kpts.tolist(),
                "energies": bstruct.energies.tolist(),
            }
            with open("tests/CrPS4/expected_band_structure.json", "w") as f:
                json.dump(data, f)

        if self.do_plot:
            emax=20.0
            bstruct.plot(emin=0.0, emax=emax)
            plt.ylabel("energy [meV]")
            plt.show()

        with open("tests/CrPS4/expected_band_structure.json", "r") as f:
            data = json.load(f)

        reference_kpts = np.array(data["kpts"])
        reference_energies = np.array(data["energies"])

        self.assertTrue(band_structures_are_equal(bstruct.path.kpts, bstruct.energies, reference_kpts, reference_energies))

    def test_primitive_cell_full_exchange_spectrum(self):
        """
        Test the band structure calculation after converting to a primitive cell. The primitive cell has
        2 bands rather than 4, but along the same k-path as the original cell the band energies should
        coincide with some bands in the larger cell.
        """

        # First calculate the non primitive spectrum
        atoms, interactions = magnon.input.create_interacting_system(
            'tests/CrPS4/CrPS4.poscar', 'tests/CrPS4/CrPS4.moments',
            'tests/CrPS4/CrPS4.exchange', 1)
        interactions = interactions.symmetrize(atoms)

        special_kpts = {
            "G": [0, 0, 0],
            "X": [1, 0, 0],
            "V1": [1, 1, 0],
            "L1": [1, 1, 1],
            "L2": [0.5, 0.5, 0.5],
            "V2": [0.5, 0.5, 0.0]
        }

        path = atoms.get_cell().bandpath(special_kpts, npoints=50, special_points=special_kpts)

        magnon_spectrum = magnon.MagnonSpectrum(atoms, interactions, num_threads=2)
        bstruct = magnon_spectrum.get_band_structure(path)

        energies=bstruct.energies.squeeze(axis=0).T

        # Now do the primitive calculation
        prim_atoms, prim_interactions, _ = magnon.build.build_primitive_cell(atoms, interactions)
        prim_magnon = magnon.MagnonSpectrum(prim_atoms, prim_interactions, num_threads=4)
        prim_bstruct = prim_magnon.get_band_structure(path)

        prim_energies=prim_bstruct.energies.squeeze(axis=0).T


        if self.do_plot:
            emax = 20.0
            prim_bstruct.plot(emin=0.0, emax=emax)
            plt.ylabel("energy [meV]")
            plt.show()

        # Check that the primitive bands exists in the non-primitive band structure. We can't compare complete bands
        # because the energies are sorted at each k-point, so we just check if, at a give k-point, the primitive
        # band energy exists in one of the non-primitive bands.
        for prim_band in prim_energies:
            self.assertTrue(
                all(np.any(np.isclose(prim_band[i], energies[:, i])) for i in range(len(prim_band)))
            )

        if self.do_overwrite:
            data = {
                "kpts": prim_bstruct.path.kpts.tolist(),
                "energies": prim_bstruct.energies.tolist(),
            }
            with open("tests/CrPS4/expected_primitive_band_structure.json", "w") as f:
                json.dump(data, f)

        if self.do_plot:
            emax = 20.0
            prim_bstruct.plot(emin=0.0, emax=emax)
            plt.ylabel("energy [meV]")
            plt.show()

        with open("tests/CrPS4/expected_primitive_band_structure.json", "r") as f:
            data = json.load(f)

        reference_kpts = np.array(data["kpts"])
        reference_energies = np.array(data["energies"])

        self.assertTrue(
            band_structures_are_equal(prim_bstruct.path.kpts, prim_bstruct.energies, reference_kpts, reference_energies))

