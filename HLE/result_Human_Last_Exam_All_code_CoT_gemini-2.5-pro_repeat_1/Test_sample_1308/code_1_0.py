import numpy as np

def calculate_glyoxal_energies():
    """
    Calculates and prints the 4 pi-electron energies of glyoxal using
    modified Hückel theory parameters h_O=1 and k_CO=0.8.
    """
    # Define the modified Hückel matrix for glyoxal (O1-C2-C3-O4)
    # H_ij = k_ij for i != j (bonded atoms)
    # H_ii = h_i
    # h_O = 1.0, h_C = 0.0
    # k_CO = 0.8, k_CC = 1.0
    huckel_matrix = np.array([
        [1.0, 0.8, 0.0, 0.0],
        [0.8, 0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0, 0.8],
        [0.0, 0.0, 0.8, 1.0]
    ])

    # The eigenvalues lambda of this matrix relate to the energy E by E = alpha + lambda*beta
    eigenvalues = np.linalg.eigvals(huckel_matrix)

    # Sort eigenvalues in descending order. Since beta is a negative energy,
    # a larger eigenvalue corresponds to a lower (more stable) energy level.
    sorted_eigenvalues = sorted(eigenvalues, reverse=True)

    print("The 4 pi-electron energies of glyoxal are:")
    # Print the energy levels in the form E = alpha +/- x*beta
    for i, lam in enumerate(sorted_eigenvalues):
        if lam < 0:
            # For negative eigenvalues, the equation is E = alpha - |lambda|*beta
            print(f"E{i+1} = α - {-lam:.3f}β")
        else:
            # For positive eigenvalues, the equation is E = alpha + lambda*beta
            print(f"E{i+1} = α + {lam:.3f}β")

# Run the calculation
calculate_glyoxal_energies()