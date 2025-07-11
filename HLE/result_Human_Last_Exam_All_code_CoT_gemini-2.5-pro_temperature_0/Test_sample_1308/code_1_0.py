import numpy as np

def solve_huckel_glyoxal():
    """
    Calculates the 4 pi-electron energy levels of glyoxal using modified Huckel theory.
    """
    # Given parameters for the heteroatom (Oxygen)
    h_O = 1.0
    k_O = 0.8

    # The Hückel matrix is constructed based on the connectivity O(1)-C(2)-C(3)-O(4)
    # A_ii = h_i
    # A_ij = k_ij for bonded atoms
    # The eigenvalues x of this matrix give the energies E = α + xβ.
    huckel_matrix = np.array([
        [h_O, k_O, 0.0, 0.0],  # Row for O(1)
        [k_O, 0.0, 1.0, 0.0],  # Row for C(2)
        [0.0, 1.0, 0.0, k_O],  # Row for C(3)
        [0.0, 0.0, k_O, h_O]   # Row for O(4)
    ])

    # Calculate the eigenvalues of the matrix
    eigenvalues = np.linalg.eigvals(huckel_matrix)

    # Sort the eigenvalues in ascending order (from most stable to least stable energy)
    eigenvalues.sort()

    print("The four pi-electron energy levels for glyoxal are:")
    # Print the energy levels in the standard format E = α + xβ
    for i, x in enumerate(eigenvalues):
        if x < 0:
            # For negative eigenvalues, the format is E = α - |x|β
            print(f"E{i+1} = α - {-x:.4f}β")
        else:
            # For positive eigenvalues, the format is E = α + xβ
            print(f"E{i+1} = α + {x:.4f}β")

# Execute the function to find and print the energies
solve_huckel_glyoxal()