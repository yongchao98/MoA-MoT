import numpy as np

# Step 1: Define the K-matrix for the Bosonic Integer Quantum Hall (BIQH) state.
# This matrix describes the topological order of the state formed by the composite fermions.
K_BIQH = np.array([[0, 1], [1, 0]])

# Step 2: Define the K-matrix for the flux attachment.
# The composite fermions are formed by attaching two fluxes to each original fermion (electron).
# This is represented by 2 times the identity matrix.
num_fluxes = 2
I = np.identity(2)
K_flux = num_fluxes * I

# Step 3: The final K-matrix is the sum of the two.
K_final = K_BIQH + K_flux

# Step 4: Print the process and the result.
print("The K-matrix for the composite fermion's state is K_BIQH:")
print(K_BIQH)
print("\nThe K-matrix for attaching two fluxes is K_flux = 2 * I:")
print(K_flux.astype(int))
print("\nThe final K-matrix is the sum: K_final = K_BIQH + K_flux")
print("\nThe element-wise calculation is:")

# Explicitly print the calculation for each element.
# This follows the request to "output each number in the final equation".
print(f"K_final[0,0] = K_BIQH[0,0] + K_flux[0,0] = {K_BIQH[0,0]} + {int(K_flux[0,0])} = {int(K_final[0,0])}")
print(f"K_final[0,1] = K_BIQH[0,1] + K_flux[0,1] = {K_BIQH[0,1]} + {int(K_flux[0,1])} = {int(K_final[0,1])}")
print(f"K_final[1,0] = K_BIQH[1,0] + K_flux[1,0] = {K_BIQH[1,0]} + {int(K_flux[1,0])} = {int(K_final[1,0])}")
print(f"K_final[1,1] = K_BIQH[1,1] + K_flux[1,1] = {K_BIQH[1,1]} + {int(K_flux[1,1])} = {int(K_final[1,1])}")

print("\nTherefore, the final K-matrix for the resulting fractional state is:")
print(K_final.astype(int))