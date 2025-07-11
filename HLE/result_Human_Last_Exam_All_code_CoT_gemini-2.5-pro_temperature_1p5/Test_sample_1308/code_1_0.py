import numpy as np

def calculate_glyoxal_energies():
    """
    Calculates the 4 pi-electron energy levels of glyoxal using modified Huckel theory.

    The problem is set up as an eigenvalue equation H C = E S C, which simplifies to
    A x = x where E = a + x*b. The matrix A is constructed using the given
    parameters for heteroatoms.
    """

    # Define the modified Huckel parameters
    # For Oxygen:
    h_O = 1.0
    k_CO = 0.8
    # For Carbon:
    h_C = 0.0
    k_CC = 1.0

    # Let's number the atoms: O(1)-C(2)-C(3)-O(4)
    # Construct the Huckel matrix A
    # A_ii = -h_i
    # A_ij = k_ij for adjacent atoms, 0 otherwise
    A = np.array([
        [-h_O, k_CO, 0,    0],
        [k_CO, -h_C, k_CC, 0],
        [0,    k_CC, -h_C, k_CO],
        [0,    0,    k_CO, -h_O]
    ])

    # Calculate the eigenvalues of the matrix A
    # The eigenvalues are the 'x' values in the energy equation E = a + x*b
    eigenvalues = np.linalg.eigvals(A)

    # Since beta is a negative quantity, a larger 'x' corresponds to a lower energy.
    # We sort the eigenvalues in descending order to list energies from most to least stable.
    eigenvalues_sorted = np.sort(eigenvalues)[::-1]

    print("The four pi-electron energies for glyoxal are:")
    # Print the final energy equations for each eigenvalue
    for i, x in enumerate(eigenvalues_sorted):
        sign = "+" if x >= 0 else "-"
        # Outputting each component of the final energy equation as requested
        print(f"E_{i+1} = α {sign} {abs(x):.4f}β")

# Run the calculation and print the results
calculate_glyoxal_energies()