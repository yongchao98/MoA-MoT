import numpy as np

# Step 1: Define the K-matrix for the nu=2 bosonic state.
K_boson = np.array([
    [0, 1],
    [1, 0]
])

# Step 2: Construct the K-matrix for the underlying unpaired composite fermions.
# This is done by the block-diagonal construction K_fermion = [[K_boson, 0], [0, -K_boson]].
K_fermion = np.block([
    [K_boson, np.zeros_like(K_boson)],
    [np.zeros_like(K_boson), -K_boson]
])

# Step 3: Define the number of attached fluxes and the J matrix.
m = 2
# The dimension of our system is now 4x4.
dim = K_fermion.shape[0]
J = np.ones((dim, dim))

# Step 4: Reverse the flux attachment transformation to find the final K-matrix.
# The formula is K_final = (K_fermion_inv - m*J)^-1.

# First, calculate the inverse of K_fermion.
# For this specific K_fermion, it is its own inverse.
K_fermion_inv = np.linalg.inv(K_fermion)

# Now, calculate the matrix inside the parentheses.
M = K_fermion_inv - m * J

# Finally, calculate the inverse of M to get the final K-matrix.
K_final = np.linalg.inv(M)

# Print the final K-matrix. The result contains integers, so we cast it for clarity.
print("The final K-matrix is:")
print(K_final.astype(int))
