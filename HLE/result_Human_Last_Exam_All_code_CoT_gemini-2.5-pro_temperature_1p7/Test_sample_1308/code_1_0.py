import numpy as np

def solve_glyoxal_huckel():
    """
    Calculates the 4 Huckel energy levels for glyoxal using modified parameters.
    """
    # Define the modified Huckel parameters
    h_O = 1.0
    k_CO = 0.8
    
    # Standard parameters for carbon
    h_C = 0.0
    k_CC = 1.0
    
    # Construct the Huckel matrix (K) for O1-C2-C3-O4
    # The matrix is defined as K_ij = (H_ij - a*delta_ij)/b
    K_matrix = np.array([
        [h_O,  k_CO,  0.0,  0.0],   # Row for O1
        [k_CO, h_C,   k_CC, 0.0],   # Row for C2
        [0.0,  k_CC,  h_C,  k_CO],  # Row for C3
        [0.0,  0.0,   k_CO, h_O]    # Row for O4
    ])

    # Calculate the eigenvalues of the matrix K. These are the 'x' values in E = a + x*b.
    eigenvalues = np.linalg.eigvals(K_matrix)
    
    # Sort the eigenvalues in descending order.
    # A larger 'x' corresponds to a lower (more stable) energy level because beta is negative.
    sorted_eigenvalues = np.sort(eigenvalues)[::-1]
    
    # Print the final energy levels in the specified format
    print("The 4 energy levels for glyoxal are found from the equation E = \u03B1 + x\u03B2.")
    print("The calculated energies are:")
    
    # Output each energy level with its corresponding eigenvalue 'x'
    # We use string formatting to insert the numbers into the equation.
    print("E1 = \u03B1 + {0:.3f}\u03B2".format(sorted_eigenvalues[0]))
    print("E2 = \u03B1 + {0:.3f}\u03B2".format(sorted_eigenvalues[1]))
    print("E3 = \u03B1 {0:+.3f}\u03B2".format(sorted_eigenvalues[2]))
    print("E4 = \u03B1 {0:+.3f}\u03B2".format(sorted_eigenvalues[3]))

# Execute the function to find and print the energies
solve_glyoxal_huckel()