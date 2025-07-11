import numpy as np

# The Hückel matrix for glyoxal O(1)-C(2)-C(3)-O(4) is constructed
# using the provided parameters: h_O = 1, k_CO = 0.8.
# For carbon, the standard parameters are h_C = 0 and k_CC = 1.
# The matrix H' has elements H'_{ii} = h_i and H'_{ij} = k_{ij}.
huckel_matrix = np.array([
    [1.0, 0.8, 0.0, 0.0],
    [0.8, 0.0, 1.0, 0.0],
    [0.0, 1.0, 0.0, 0.8],
    [0.0, 0.0, 0.8, 1.0]
])

# The eigenvalues (λ) of this matrix determine the energy levels.
eigenvalues = np.linalg.eigvals(huckel_matrix)

# The energy levels are E = α + λβ. Since β is a negative value,
# a larger eigenvalue corresponds to a lower (more stable) energy.
# We sort the eigenvalues in descending order to list energies from lowest to highest.
sorted_eigenvalues = np.sort(eigenvalues)[::-1]

# Print the final energy equations.
print("The four Hückel energies for glyoxal are:")
for i, lam in enumerate(sorted_eigenvalues):
    sign = "+" if lam >= 0 else "-"
    abs_lam = abs(lam)
    # The final equation includes all components: α, the sign, the calculated value, and β.
    print(f"E{i+1} = α {sign} {abs_lam:.4f}β")