import numpy as np

def calculate_glyoxal_energies():
    """
    Calculates the 4 pi-electron energies of glyoxal using modified Huckel theory.
    """
    # The Huckel matrix H is constructed for the O(1)-C(2)-C(3)-O(4) system.
    # The matrix elements are determined by the parameters h and k:
    # H_ii = h_i
    # H_ij = k_ij (for bonded i,j)
    # Parameters: h_O = 1, h_C = 0, k_CO = 0.8, k_CC = 1.
    H = np.array([
        [1.0, 0.8, 0.0, 0.0],
        [0.8, 0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0, 0.8],
        [0.0, 0.0, 0.8, 1.0]
    ])

    # The eigenvalues of this matrix give the energy coefficients x in E = a + x*b.
    eigenvalues = np.linalg.eigvals(H)

    # Sort the eigenvalues to present the energies from lowest to highest.
    eigenvalues.sort()

    print("The 4 pi-electron energy levels for glyoxal are:")
    # Print the energy levels in terms of alpha (a) and beta (b).
    # The format E = a + (x)b is used for clarity.
    for i, eigval in enumerate(eigenvalues):
        print(f"E{i+1} = α + ({eigval:.4f})β")

# Execute the function to print the results
calculate_glyoxal_energies()