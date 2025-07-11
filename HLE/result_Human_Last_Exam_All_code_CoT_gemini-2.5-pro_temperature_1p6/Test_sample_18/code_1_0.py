import numpy as np

# This script calculates the K-matrix for a fractional quantum Hall state
# based on the properties of a related bosonic state.

# Plan:
# 1. Define the initial K-matrix for the bosons (K_B), which is the Pauli sigma_x matrix.
# 2. The physical description implies that the bosons are Cooper pairs of composite fermions (CF).
#    To find the K-matrix for the CFs (K_CF), we must reverse the pairing process.
# 3. This "un-pairing" involves two steps:
#    a. Charge transformation: A boson's charge is double a CF's charge. This results in
#       transforming the K-matrix by a factor of 4 (K' = 4 * K_B).
#    b. Statistics transformation: To change the particles from bosons to fermions, we add
#       the identity matrix I (K_CF = K' + I).
# 4. The final formula is K_CF = 4 * K_B + I.
# 5. The script will compute and print the result, showing each number in the equation.

# Step 1: Define the initial K-matrix for the bosons
K_B = np.array([[0, 1],
                [1, 0]])

# Step 2: Define the identity matrix for the statistics transformation
I = np.identity(2, dtype=int)

# Step 3: Apply the transformation K_CF = 4 * K_B + I
K_CF = 4 * K_B + I

# Step 4: Print the detailed calculation as requested
print("The K-matrix for the initial bosonic state is K_B = sigma_x:")
print(K_B)
print("\nThe K-matrix for the composite fermion state, K_CF, is derived from K_B.")
print("The calculation is K_CF = 4 * K_B + I.")
print("\nHere is the step-by-step calculation with all numbers shown:")

# Print the equation with the initial matrix values
print(f"K_CF = 4 * [[{K_B[0, 0]}, {K_B[0, 1]}], [{K_B[1, 0]}, {K_B[1, 1]}]] + [[{I[0, 0]}, {I[0, 1]}], [{I[1, 0]}, {I[1, 1]}]]")

# Print the result of the multiplication
K_prime = 4 * K_B
print(f"     = [[{K_prime[0, 0]}, {K_prime[0, 1]}], [{K_prime[1, 0]}, {K_prime[1, 1]}]] + [[{I[0, 0]}, {I[0, 1]}], [{I[1, 0]}, {I[1, 1]}]]")

# Print the final resulting matrix
print(f"     = [[{K_CF[0, 0]}, {K_CF[0, 1]}], [{K_CF[1, 0]}, {K_CF[1, 1]}]]")

print("\nThus, the final K-matrix is:")
print(K_CF)