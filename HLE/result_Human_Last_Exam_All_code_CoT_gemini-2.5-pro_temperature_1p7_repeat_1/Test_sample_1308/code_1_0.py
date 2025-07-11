import numpy as np

# Define the Hückel theory parameters for the heteroatoms
# hO is for the oxygen coulomb integral: α_O = α + h_O * β
# kO is for the carbon-oxygen resonance integral: β_CO = k_O * β
hO = 1.0
kO = 0.8

# Construct the Hückel parameter matrix M' for glyoxal (O-C-C-O).
# The energy levels are given by E = α + λβ, where λ are the eigenvalues of this matrix.
# The structure is O(1)-C(2)-C(3)-O(4)
M_prime = np.array([
    [hO,  kO,  0.0, 0.0],
    [kO,  0.0, 1.0, 0.0],
    [0.0, 1.0, 0.0, kO ],
    [0.0, 0.0, kO,  hO ]
])

# Calculate the eigenvalues of the matrix
eigenvalues = np.linalg.eigvals(M_prime)

# Since β is a negative quantity, a larger (more positive) eigenvalue
# corresponds to a lower (more stable) energy level.
# We sort the eigenvalues in descending order to list energies from lowest to highest.
sorted_eigenvalues = np.sort(eigenvalues)[::-1]

# Print the final energy levels in the standard Huckel theory format.
# The code explicitly shows the numerical coefficient for β in each equation.
print("The 4 pi-electron energy levels (E) of glyoxal are:")
for i, lam in enumerate(sorted_eigenvalues):
    if lam < 0:
        # For negative lambda, print as 'α - |λ|β' for clarity
        print(f"E{i+1} = α - {-lam:.4f}β")
    else:
        # For positive lambda, print as 'α + λβ'
        print(f"E{i+1} = α + {lam:.4f}β")