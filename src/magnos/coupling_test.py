# import numpy.testing as npt
#
# from unittest import TestCase
#
# from .interactions import InteractionList, apply_bond_reversal_symmetry
#
#
# class IntrinsicSymmetryTests(TestCase):
#
#     def test_intrinsic_symmetry_fails_with_inconsistent_data(self):
#         """
#         Test that an exception is raised when the length of the input data is inconsistent.
#         :return:
#         """
#
#         # Length of the input is incorrect
#
#         coupling_indices = [[0, 0], [0, 0]]
#         coupling_vectors = [[0.0, 0.0, 1.0]]
#         coupling_matrices = [[[1.0, 2.0, 3.0], [-2.0, 1.0, 4.0], [-3.0, -4.0, 1.0]]]
#
#         with self.assertRaises(Exception):
#             apply_bond_reversal_symmetry(coupling_indices, coupling_vectors, coupling_matrices)
#
#         coupling_indices = [[0, 0]]
#         coupling_vectors = [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]
#         coupling_matrices = [[[1.0, 2.0, 3.0], [-2.0, 1.0, 4.0], [-3.0, -4.0, 1.0]]]
#
#         with self.assertRaises(Exception):
#             apply_bond_reversal_symmetry(coupling_indices, coupling_vectors, coupling_matrices)
#
#         coupling_indices = [[0, 0]]
#         coupling_vectors = [[0.0, 0.0, 1.0]]
#         coupling_matrices = [[[1.0, 2.0, 3.0], [-2.0, 1.0, 4.0], [-3.0, -4.0, 1.0]], [[1.0, 2.0, 3.0], [-2.0, 1.0, 4.0], [-3.0, -4.0, 1.0]]]
#
#         with self.assertRaises(Exception):
#             apply_bond_reversal_symmetry(coupling_indices, coupling_vectors, coupling_matrices)
#
#         # Multiply defined interactions
#
#         coupling_indices = [[0, 0], [0, 0], [0, 0]]
#         coupling_vectors = [[0.0, 0.0, 1.0], [0.0, 0.0,-1.0], [0.0, 0.0,-1.0]]
#         coupling_matrices = [
#             [[1.0, 2.0, 3.0], [-2.0, 1.0, 4.0], [-3.0, -4.0, 1.0]],
#             [[1.0, 2.0, 3.0], [-2.0, 1.0, 4.0], [-3.0, -4.0, 1.0]],
#             [[1.0, 2.0, 3.0], [-2.0, 1.0, 4.0], [-3.0, -4.0, 1.0]]]
#
#         with self.assertRaises(RuntimeError):
#             apply_bond_reversal_symmetry(coupling_indices, coupling_vectors, coupling_matrices)
#
#
#     def test_intrinsic_symmetry_generates_new_couplings(self):
#         """
#         Test the correct new indices, vectors and matrices are generated
#         """
#
#
#         coupling_indices = [[0, 0]]
#         coupling_vectors = [[0.0, 0.0, 1.0]]
#         coupling_matrices = [[[1.0, 2.0, 3.0], [-2.0, 1.0, 4.0], [-3.0, -4.0, 1.0]]]
#
#         new_indices, new_vectors, new_matrices = apply_bond_reversal_symmetry(coupling_indices,
#                                                                               coupling_vectors,
#                                                                               coupling_matrices)
#
#         npt.assert_array_equal(new_indices, [[0, 0]])
#         npt.assert_array_equal(new_vectors, [[0.0, 0.0, -1.0]])
#         npt.assert_array_equal(new_matrices, [[[1.0, -2.0, -3.0], [2.0, 1.0, -4.0], [3.0, 4.0, 1.0]]])
#
#
#         coupling_indices = [[0, 1]]
#         coupling_vectors = [[0.0, 0.0, 1.0]]
#         coupling_matrices = [[[1.0, 2.0, 3.0], [-2.0, 1.0, 4.0], [-3.0, -4.0, 1.0]]]
#
#         new_indices, new_vectors, new_matrices = apply_bond_reversal_symmetry(coupling_indices,
#                                                                               coupling_vectors,
#                                                                               coupling_matrices)
#
#         npt.assert_array_equal(new_indices, [[1, 0]])
#         npt.assert_array_equal(new_vectors, [[0.0, 0.0, -1.0]])
#         npt.assert_array_equal(new_matrices, [[[1.0, -2.0, -3.0], [2.0, 1.0, -4.0], [3.0, 4.0, 1.0]]])
#
#
#     def test_intrinsic_symmetry_fails_with_broken_coupling_symmetry(self):
#         """
#         Test that an exception is raised when the intrinsic symmetry is broken by J_ij != J_ji^T
#         """
#         with self.assertRaises(RuntimeError):
#
#             coupling_indices = [[0, 1], [1, 0]]
#             coupling_vectors = [[0.0, 0.0, 1.0], [0.0, 0.0,-1.0]]
#             coupling_matrices = [[[1.0, 2.0, 3.0], [-2.0, 1.0, 4.0], [-3.0, -4.0, 1.0]],
#                                  [[1.0, 2.0, 3.0], [-2.0, 1.0, 4.0], [-3.0, -4.0, 1.0]]]
#
#             new_indices, new_vectors, new_matrices = apply_bond_reversal_symmetry(coupling_indices,
#                                                                                   coupling_vectors,
#                                                                                   coupling_matrices)
#
