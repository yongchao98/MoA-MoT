import numpy as np

# Define the matrices a and b
a = np.array([[-21, 242], [-2, 23]])
b = np.array([[-19, 200], [-2, 21]])

# Define the generators of Gamma(2)
A = np.array([[1, 2], [0, 1]])
B = np.array([[1, 0], [2, 1]])

# Inverse of A and B
A_inv = np.linalg.inv(A).astype(int)
B_inv = np.linalg.inv(B).astype(int)

# --- Decomposition of a ---
# M1 = A^{-5} * a
A_inv_5 = np.linalg.matrix_power(A_inv, 5)
M1_a = A_inv_5 @ a
# M2 = B^{-1} * M1
M2_a = B_inv @ M1_a

# Check M2_a structure
Z = np.array([[-1, 0], [0, -1]])
A_inv_6 = np.linalg.matrix_power(A_inv, 6)
# Z * A^{-6}
check_M2_a = Z @ A_inv_6

# The abelianization vectors
v_a = np.array([5, 0]) + np.array([0, 1]) - np.array([6, 0])

# --- Decomposition of b ---
# N1 = A^{-5} * b
N1_b = A_inv_5 @ b
# N2 = B * N1
N2_b = B @ N1_b

# Check N2_b structure
A_inv_5_check = np.linalg.matrix_power(A_inv, 5)

# The abelianization vectors
v_b = np.array([5, 0]) - np.array([0, 1]) - np.array([5, 0])

# --- Calculate the index ---
# Index of Gamma(2) in SL_2(Z)
index_G_Gamma2 = 6

# Matrix of abelianization vectors
M_ab = np.column_stack((v_a, v_b))

# Index of H in Gamma(2)
index_Gamma2_H = int(abs(np.linalg.det(M_ab)))

# Total index [G:H]
total_index = index_G_Gamma2 * index_Gamma2_H

# Print the results
print("Matrix a:")
print(a)
print("\nMatrix b:")
print(b)
print("\nAbelianization of a: (exp_A, exp_B) = ({}, {})".format(v_a[0], v_a[1]))
print("Abelianization of b: (exp_A, exp_B) = ({}, {})".format(v_b[0], v_b[1]))
print("\nMatrix of abelianized vectors:")
print(M_ab)
print("\nIndex [Gamma(2) : H] = |det(M_ab)| = {}".format(index_Gamma2_H))
print("Index [G : Gamma(2)] = {}".format(index_G_Gamma2))
print("\nFinal Index [G : H] = [G : Gamma(2)] * [Gamma(2) : H] = {} * {} = {}".format(index_G_Gamma2, index_Gamma2_H, total_index))
print(f"Final answer: {total_index}")
