import numpy as np

# Step 1: Define the K-matrix for the BIQH state.
K_boson = np.array([[0, 1], 
                    [1, 0]])

# Step 2: Define the transformation matrix W for Cooper pairing.
# l_1 = CF_1 + CF_2  -> column vector [1, 1]
# l_2 = CF_1 - CF_2  -> column vector [1, -1]
W = np.array([[1, 1], 
              [1, -1]])

# Step 3: Calculate the K-matrix for the composite fermions (K_CF).
# The transformation is K_CF = W^T * K_boson * W
WT = W.transpose()
K_CF = WT @ K_boson @ W

# Step 4: Define the K-matrix for the attached fluxes.
# "Two fluxes" corresponds to 2*I.
K_flux = 2 * np.identity(2)

# Step 5: Calculate the final K-matrix for the fermionic state.
# K_final = K_CF - K_flux
K_final = K_CF - K_flux

# Print the results of the calculation step-by-step
print("The K-matrix for the nu=2 BIQH state is:")
print(f"K_boson = \n{K_boson}\n")

print("The transformation matrix for Cooper pairing is:")
print(f"W = \n{W}\n")

print("The K-matrix for the composite fermions is K_CF = W^T * K_boson * W:")
print(f"K_CF = \n{K_CF}\n")

print("The final K-matrix for the resulting fractional state is K_final = K_CF - 2*I:")
print("K_final = [[{:.0f}, {:.0f}], [{:.0f}, {:.0f}]] - [[{:.0f}, {:.0f}], [{:.0f}, {:.0f}]] = [[{:.0f}, {:.0f}], [{:.0f}, {:.0f}]]".format(
    K_CF[0, 0], K_CF[0, 1], K_CF[1, 0], K_CF[1, 1],
    K_flux[0, 0], K_flux[0, 1], K_flux[1, 0], K_flux[1, 1],
    K_final[0, 0], K_final[0, 1], K_final[1, 0], K_final[1, 1]
))
