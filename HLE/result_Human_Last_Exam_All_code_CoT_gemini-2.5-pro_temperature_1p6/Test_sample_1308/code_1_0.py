import numpy as np

def calculate_glyoxal_energies():
    """
    Calculates the 4 pi system energies of glyoxal using modified Huckel theory.
    """
    # Define the Huckel theory parameters for oxygen
    h_O = 1.0
    k_CO = 0.8

    # Define the parameters for carbon (standard)
    h_C = 0.0
    k_CC = 1.0

    # Construct the Huckel matrix for glyoxal (O1-C2-C3-O4)
    # The matrix elements are H'_ii = h_x and H'_ij = k_ij
    H_prime = np.array([
        [h_O,  k_CO, 0.0,  0.0],  # Atom 1 (O)
        [k_CO, h_C,  k_CC, 0.0],  # Atom 2 (C)
        [0.0,  k_CC, h_C,  k_CO], # Atom 3 (C)
        [0.0,  0.0,  k_CO, h_O]   # Atom 4 (O)
    ])

    # The eigenvalues 'x' of this matrix relate to the energy levels 'E'
    # by the equation E = α + xβ.
    # We use np.linalg.eigvalsh as the matrix is symmetric.
    eigenvalues = np.linalg.eigvalsh(H_prime)

    # The resonance integral β is negative. Therefore, to get the energies
    # from lowest to highest, we must sort the eigenvalues from largest to smallest.
    sorted_eigenvalues = np.sort(eigenvalues)[::-1]

    print("The four π molecular orbital energies for glyoxal are:")

    # Print each energy level in the standard format
    for i, x in enumerate(sorted_eigenvalues):
        sign = '+' if x >= 0 else '-'
        # The prompt requires outputting the numbers in the final equation
        # Example: E = α + (1.882)β or E = α - (1.229)β
        # The code below prints it as E = α + 1.882β or E = α + -1.229β which is also clear.
        print(f"E{i+1} = α + ({x:.3f})β")

# Execute the function
calculate_glyoxal_energies()