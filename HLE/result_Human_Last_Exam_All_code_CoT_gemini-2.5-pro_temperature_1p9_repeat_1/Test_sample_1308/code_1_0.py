import numpy as np

def find_glyoxal_energies():
    """
    Calculates the 4 pi-electron energies of glyoxal using modified Huckel theory.

    The secular equation is set up for the O(1)-C(2)-C(3)-O(4) pi system.
    The energy levels E are given by E = α + xβ, where x represents the
    eigenvalues of the Huckel matrix.

    Given parameters for the heteroatom (Oxygen):
    h_O = 1
    k_O = 0.8
    """

    # Define the modified Huckel theory parameters for Oxygen
    h_O = 1.0
    k_CO = 0.8

    # Construct the Huckel matrix for glyoxal: O(1)-C(2)-C(3)-O(4)
    # The matrix elements are determined by the connectivity and atom types.
    # A[i,i] = h_X for a heteroatom, 0 for carbon.
    # A[i,j] = k_XY for a bond, 1 for C-C bond, 0 if not bonded.
    huckel_matrix = np.array([
        [h_O,  k_CO, 0.0,  0.0],  # Atom 1 (O) connected to Atom 2 (C)
        [k_CO, 0.0,  1.0,  0.0],  # Atom 2 (C) connected to Atom 1 (O) and 3 (C)
        [0.0,  1.0,  0.0,  k_CO], # Atom 3 (C) connected to Atom 2 (C) and 4 (O)
        [0.0,  0.0,  k_CO, h_O]   # Atom 4 (O) connected to Atom 3 (C)
    ])

    # The solutions 'x' to the secular determinant are the eigenvalues of this matrix.
    eigenvalues = np.linalg.eigvals(huckel_matrix)

    # Sort the eigenvalues in ascending order to list energies from lowest to highest.
    eigenvalues.sort()

    print("The four pi-electron energies of glyoxal are:")

    # Print the energy levels. The energy E is expressed as E = α + xβ.
    # Each number in the final equation (the coefficient of β) is printed.
    for i, x in enumerate(eigenvalues):
        if x < 0:
            # For negative eigenvalues, format the equation as E = α - |x|β
            print(f"E{i+1} = α - {-x:.4f}β")
        else:
            # For positive eigenvalues, format the equation as E = α + xβ
            print(f"E{i+1} = α + {x:.4f}β")

# Run the calculation and print the results
find_glyoxal_energies()