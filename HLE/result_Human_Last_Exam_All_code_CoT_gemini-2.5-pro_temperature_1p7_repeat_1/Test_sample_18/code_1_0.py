import numpy as np

# Step 1: Define the K-matrix for the bosonic state.
# The K-matrix for the nu=2 bosonic IQH state is given as the Pauli sigma_x matrix.
K_b = np.array([[0, 1],
                [1, 0]])

# Step 2: Determine the K-matrix of the Composite Fermions (CFs).
# The bosons are s-wave Cooper pairs of CFs. The transformation is K_b = 4 * K_CF.
# To find the K-matrix of the CFs, we reverse this transformation.
K_CF = K_b / 4.0

# Step 3: Determine the K-matrix of the final electronic state.
# The CFs are formed by attaching m=2 fluxes to each fundamental fermion.
# The transformation is K_e = K_CF + m * I, where I is the identity matrix.
m = 2
I = np.identity(2)
K_final = K_CF + m * I

# Step 4: Print the final K-matrix.
print("The K-matrix of the resulting fractional state is:")
print(f"K = [[{K_final[0, 0]}, {K_final[0, 1]}],")
print(f"     [{K_final[1, 0]}, {K_final[1, 1]}]]")
