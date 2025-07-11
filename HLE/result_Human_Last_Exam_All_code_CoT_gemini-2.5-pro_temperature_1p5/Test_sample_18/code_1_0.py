import numpy as np

# Step 1: Define the initial K-matrix for the bosonic state.
K_b = np.array([[0, 1], 
                [1, 0]])

# Step 2: Calculate the K-matrix for the composite fermions (K_cf) by "un-pairing" the bosons.
# The transformation is K_cf = 4 * K_b.
K_cf = 4 * K_b

# Step 3: Define the flux attachment matrix.
# Two fluxes are attached to each of the two types of fermions.
# This corresponds to 2 * I, where I is the identity matrix.
K_flux = 2 * np.identity(2)

# Step 4: Calculate the final K-matrix for the fundamental fermions (K_f).
# The relation is K_cf = K_f + K_flux, so K_f = K_cf - K_flux.
K_f = K_cf - K_flux

# Step 5: Print the results, including the final equation with all numbers.
print("The K-matrix for the initial bosonic state (K_b) is:")
print(K_b)
print("\nThe K-matrix for the constituent composite fermions (K_cf = 4 * K_b) is:")
print(K_cf)
print("\nThe final K-matrix for the fundamental fermions (K_f) is found by solving K_f = K_cf - K_flux.")
print("\nThe final equation is:")
# Use .tolist() for a cleaner print format of the equation
print(f"{K_cf.astype(int).tolist()} - {K_flux.astype(int).tolist()} = {K_f.astype(int).tolist()}")

print("\nTherefore, the resulting K-matrix is:")
print(K_f.astype(int))