# import unittest
# import numpy as np
#
# from magnon.build import cell_transformation, convert_kpts
# from magnon.lattice import cartesian_k_from_scaled, reciprocal_lattice, scaled_k_from_cartesian
#
# class TestCellTransformation(unittest.TestCase):
#
#     def test_simple_transform(self):
#         original_cell = np.array(
#             [[1, 0, 0],
#             [0, 1, 0],
#             [0, 0, 1]])
#
#         theta = np.pi / 4  # 45 degrees in radians
#         Rz = np.array([
#             [np.cos(theta), -np.sin(theta), 0],
#             [np.sin(theta), np.cos(theta), 0],
#             [0, 0, 1]])
#
#         transformed_cell = Rz @ original_cell
#         transformation = cell_transformation(original_cell, transformed_cell)
#
#         self.assertTrue(np.allclose(Rz, transformation))
#
# class TestKSpaceConversion(unittest.TestCase):
#
#     def test_kspace_round_trip(self):
#         # Use a monoclinic system to avoid accidentally passing due to high symmetry.
#         # This is from the CrPS4 example, lattice vectors are the rows of the matrix
#         # which is the convention used in ASE and spglib (python API!).
#         cell = np.array([
#             [   1.0, 0.0,    0.0],
#             [   0.0, 0.6687, 0.0],
#             [ -0.02, 0.0,    1.1246]])
#
#         recip_lattice = reciprocal_lattice(cell)
#
#         # Some high symmetry points
#         k_scaled = [
#             [0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1], [0.5, 0.5, 0.5], [0.5, 0.5, 0.0]]
#
#         k_cartesian = cartesian_k_from_scaled(k_scaled, recip_lattice)
#         k_scaled2 = scaled_k_from_cartesian(k_cartesian, recip_lattice)
#
#         self.assertTrue(np.allclose(k_scaled, k_scaled2))
#
#     def test_simple_kpath_conversion(self):
#         """
#         Rotate a simple cubic cell by 45 degrees around the z-axis. Check that converting the high symmetry points in
#         the original reciprocal cell gives the correct points in the new reciprocal cell.
#         Returns:
#
#         """
#
#         original_cell = np.array(
#             [[1, 0, 0],
#             [0, 1, 0],
#             [0, 0, 1]])
#
#         theta = np.pi / 4  # 45 degrees in radians
#         Rz = np.array([
#             [np.cos(theta), -np.sin(theta), 0],
#             [np.sin(theta), np.cos(theta), 0],
#             [0, 0, 1]])
#
#         transformed_cell = Rz @ original_cell
#
#         kpts = {
#             "G": [0, 0, 0],
#             "M": [0.5, 0.5, 0],
#             "R": [0.5, 0.5, 0.5],
#             "X": [0.0, 0.5, 0.0]
#         }
#
#         converted_kpts = {
#             symbol: convert_kpts(np.array(kpoint), original_cell, transformed_cell)
#             for symbol, kpoint in kpts.items()
#         }
#
#         self.assertTrue(np.allclose(converted_kpts['G'], np.array([0.0, 0.0, 0])))
#         self.assertTrue(np.allclose(converted_kpts['M'], np.array([1.0/np.sqrt(2), 0, 0])))
#         self.assertTrue(np.allclose(converted_kpts['R'], np.array([1.0 / np.sqrt(2), 0, 0.5])))
#         self.assertTrue(np.allclose(converted_kpts['X'], np.array([np.sqrt(1/8), np.sqrt(1/8), 0.0])))
#
#     def test_volume_change_kpath_conversion(self):
#         """
#         Rotate a simple cubic cell by 45 degrees around the z-axis. Check that converting the high symmetry points in
#         the original reciprocal cell gives the correct points in the new reciprocal cell.
#
#         Original lattice vectors:
#         a1 = a(1, 0, 0)
#         a2 = a(0, 1, 0)
#         a3 = a(0, 0, 1)
#
#         Original reciprocal vectors:
#         b1 = (2π/a)(1, 0, 0)
#         b2 = (2π/a)(0, 1, 0)
#         b3 = (2π/a)(0, 0, 1)
#
#         If we rotate the lattice by 45 degrees around the z-axis we get
#
#         New lattice vectors:
#         a1 = a(1/√2, -1/√2, 0)
#         a2 = a(1/√2, 1/√2, 0)
#         a3 = a(0, 0, 1)
#
#         New reciprocal vectors:
#         b1 = (2π/a)(1/√2, -1/√2, 0)
#         b2 = (2π/a)(0, 1, 0)
#         b3 = (2π/a)(0, 0, 1)
#
#         """
#
#         original_cell = np.array(
#             [[1, 0, 0],
#              [0, 1, 0],
#              [0, 0, 1]])
#
#         theta = np.pi / 4  # 45 degrees in radians
#         Rz = np.array([
#             [np.cos(theta), -np.sin(theta), 0],
#             [np.sin(theta), np.cos(theta), 0],
#             [0, 0, 1]])
#
#         transformed_cell = Rz @ original_cell
#
#         kpts = {
#             "G": [0, 0, 0],
#             "M": [0.5, 0.5, 0],
#             "R": [0.5, 0.5, 0.5],
#             "X": [0.0, 0.5, 0.0]
#         }
#
#         converted_kpts = {
#             symbol: convert_kpts(np.array(kpoint), original_cell, transformed_cell)
#             for symbol, kpoint in kpts.items()
#         }
#
#         self.assertTrue(np.allclose(converted_kpts['G'], np.array([0.0, 0.0, 0])))
#         self.assertTrue(np.allclose(converted_kpts['M'], np.array([1.0 / np.sqrt(2), 0, 0])))
#         self.assertTrue(np.allclose(converted_kpts['R'], np.array([1.0 / np.sqrt(2), 0, 0.5])))
#         self.assertTrue(np.allclose(converted_kpts['X'], np.array([np.sqrt(1 / 8), np.sqrt(1 / 8), 0.0])))