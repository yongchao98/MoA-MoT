import numpy as np
import sys

def solve_huckel_energies_for_glyoxal():
    """
    Calculates and prints the Hückel pi-electron energies for glyoxal.

    The molecule glyoxal, O=CH-CH=O, has a linear pi system O(1)-C(2)-C(3)-O(4).
    The modified Hückel parameters are:
    - h_O = 1 (for oxygen atoms 1 and 4)
    - h_C = 0 (for carbon atoms 2 and 3)
    - k_CO = 0.8 (for C-O bonds 1-2 and 3-4)
    - k_CC = 1 (for C-C bond 2-3)

    The energy E is related to the eigenvalues 'x' of the Hückel matrix by E = α - xβ.
    The 'x' values are the eigenvalues of the matrix A, where A_ii = -h_i and A_ij = -k_ij.
    """

    # Construct the Hückel matrix A
    # A = [ [-h1, -k12, -k13, -k14],
    #       [-k21, -h2, -k23, -k24],
    #       ...                       ]
    huckel_matrix = np.array([
        [-1, -0.8, 0, 0],    # Row for O(1)
        [-0.8, 0, -1, 0],    # Row for C(2)
        [0, -1, 0, -0.8],    # Row for C(3)
        [0, 0, -0.8, -1]     # Row for O(4)
    ])

    # The eigenvalues of this matrix are the 'x' values
    try:
        x_values = np.linalg.eigvals(huckel_matrix)
    except np.linalg.LinAlgError:
        print("Error: Eigenvalue computation failed.", file=sys.stderr)
        return

    # Sort the eigenvalues in descending order to list energies
    # from lowest to highest.
    x_values_sorted = np.sort(x_values)[::-1]

    print("The four pi-electron energies (E) for glyoxal are found by solving the Hückel secular determinant.")
    print("The solutions, x, relate to the energy by the equation: E = α - xβ.")
    
    print("\nThe four calculated energy levels are:")
    
    # Print each energy value, including the numbers in the final equation.
    for i, x in enumerate(x_values_sorted):
        # Format the energy expression as E = α +/- Cβ
        sign = "+"
        coefficient = -x
        if coefficient < 0:
            sign = "-"
            coefficient = abs(coefficient)
        
        print(f"E{i+1} = α {sign} {coefficient:.4f}β")

if __name__ == '__main__':
    solve_huckel_energies_for_glyoxal()
