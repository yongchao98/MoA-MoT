import numpy as np

# Step 1: Define the K-matrix for the bosonic integer quantum Hall state.
K_b = np.array([[0, 1], 
                [1, 0]])

# Step 2: The bosons are Cooper pairs of composite fermions.
# The K-matrix for the constituents (K_CF) is related to the K-matrix for the pairs (K_b)
# by K_CF = 8 * K_b.
K_CF = 8 * K_b

# Step 3: The composite fermions are formed by attaching two flux quanta to each electron.
# This corresponds to adding a Chern-Simons matrix C_2.
# C_2 is a matrix where all entries are 2.
C_2 = np.array([[2, 2], 
                [2, 2]])

# The relationship is K_CF = K_e + C_2, where K_e is the electronic K-matrix.

# Step 4: Solve for K_e.
# K_e = K_CF - C_2
K_e = K_CF - C_2

# Print the resulting K-matrix for the fractional state.
print("The K-matrix of the resulting fractional state is:")
# The final output needs to show the full equation.
# Final Equation: K_e = 8 * K_b - C_2
# K_e = 8 * [[0, 1], [1, 0]] - [[2, 2], [2, 2]]
# K_e = [[0, 8], [8, 0]] - [[2, 2], [2, 2]]
# K_e = [[-2, 6], [6, -2]]

print(f"K = 8 * [[{K_b[0,0]}, {K_b[0,1]}], [{K_b[1,0]}, {K_b[1,1]}]] - [[{C_2[0,0]}, {C_2[0,1]}], [{C_2[1,0]}, {C_2[1,1]}]]")
print(f"  = [[{K_CF[0,0]}, {K_CF[0,1]}], [{K_CF[1,0]}, {K_CF[1,1]}]] - [[{C_2[0,0]}, {C_2[0,1]}], [{C_2[1,0]}, {C_2[1,1]}]]")
print(f"  = [[{K_e[0,0]}, {K_e[0,1]}], [{K_e[1,0]}, {K_e[1,1]}]]")