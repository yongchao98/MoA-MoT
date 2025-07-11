import numpy as np

# Step 1: Define the K-matrix for the composite fermion (CF) state.
# This is the nu=2 integer quantum Hall state for a two-component system,
# which corresponds to the identity matrix.
K_CF = np.array([[1, 0],
                 [0, 1]])

# Step 2: Define the K-matrix representing the attached flux.
# m=2 fluxes are attached to each fermion, coupling to the total density.
# The coupling matrix C has all entries as 1. K_flux = m * C.
m = 2
C = np.array([[1, 1],
              [1, 1]])
K_flux = m * C

# Step 3: Calculate the final K-matrix for the electrons.
K_electron = K_CF + K_flux

# Step 4: Print the final equation clearly showing each component.
print("The K-matrix of the resulting fractional state is calculated as follows:")
print("K_electron = K_CF + K_flux\n")

print(f"  ( {K_electron[0, 0]}  {K_electron[0, 1]} )   =   ( {K_CF[0, 0]}  {K_CF[0, 1]} )   +   ( {K_flux[0, 0]}  {K_flux[0, 1]} )")
print(f"  ( {K_electron[1, 0]}  {K_electron[1, 1]} )       ( {K_CF[1, 0]}  {K_CF[1, 1]} )       ( {K_flux[1, 0]}  {K_flux[1, 1]} )\n")

print("Final K-matrix:")
print(K_electron)