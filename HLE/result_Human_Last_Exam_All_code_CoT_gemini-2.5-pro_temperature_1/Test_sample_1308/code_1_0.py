import numpy as np

def find_glyoxal_energies():
    """
    Calculates the 4 Hückel molecular orbital energies for glyoxal using
    the modified parameters hO=1 and kO=0.8.
    """
    # Define the modified Hückel parameters for the heteroatom (Oxygen)
    h_O = 1.0
    k_CO = 0.8

    # Standard parameters for Carbon atoms and C-C bonds
    h_C = 0.0
    k_CC = 1.0

    # Construct the Hückel matrix for the O(1)-C(2)-C(3)-O(4) pi system
    # The energy levels are E = α + λβ, where λ are the eigenvalues of this matrix.
    Huckel_matrix = np.array([
        [h_O,  k_CO,   0.0,    0.0 ],
        [k_CO, h_C,    k_CC,   0.0 ],
        [0.0,  k_CC,   h_C,    k_CO],
        [0.0,  0.0,    k_CO,   h_O ]
    ])

    # Calculate the eigenvalues of the matrix
    eigenvalues = np.linalg.eigvals(Huckel_matrix)

    # The energy levels are given by E = α + λβ. Since β is a negative
    # energy unit, a larger positive value of λ corresponds to a lower
    # (more stable) energy level. We sort the eigenvalues in descending
    # order to list the energies from the lowest (E1) to the highest (E4).
    sorted_eigenvalues = np.sort(eigenvalues)[::-1]

    print("The 4 Hückel molecular orbital energies for glyoxal are:")
    # Print each energy level equation, showing all the components.
    for i, val in enumerate(sorted_eigenvalues):
        # Use a different sign for positive and negative eigenvalues for clarity
        if val >= 0:
            sign = "+"
        else:
            sign = "-"
        
        # Output the final equation for each energy level
        print(f"E{i+1} = α {sign} {abs(val):.3f}β")

# Run the calculation and print the results
find_glyoxal_energies()