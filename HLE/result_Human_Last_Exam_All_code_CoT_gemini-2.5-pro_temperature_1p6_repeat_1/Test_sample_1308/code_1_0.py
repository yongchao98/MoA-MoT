import numpy as np

def calculate_glyoxal_energies():
    """
    Calculates the 4 pi-electron energy levels of glyoxal using
    modified Hückel theory.
    """
    # Parameters for heteroatoms (Oxygen)
    hO = 1.0
    kCO = 0.8

    # Standard parameters for Carbon
    hC = 0.0
    kCC = 1.0

    # Construct the Hückel matrix for glyoxal (O1-C2-C3-O4).
    # We solve for the eigenvalues 'x' in the energy equation E = α + xβ
    # by setting α=0 and β=1 in the matrix.
    huckel_matrix = np.array([
        [hO,   kCO, 0,     0  ],  # Row for O1
        [kCO,  hC,  kCC,   0  ],  # Row for C2
        [0,    kCC, hC,    kCO],  # Row for C3
        [0,    0,   kCO,   hO ]   # Row for O4
    ])

    # Calculate the eigenvalues of the Hückel matrix
    eigenvalues = np.linalg.eigvals(huckel_matrix)

    # Sort the eigenvalues in descending order, as is conventional for
    # Hückel energies (from most bonding to most anti-bonding).
    eigenvalues_sorted = np.sort(eigenvalues)[::-1]

    print("The 4 pi-electron energy levels for glyoxal are:")

    # Print the energy levels in the standard E = α ± xβ format.
    for i, x in enumerate(eigenvalues_sorted):
        # Determine the operator sign for the output equation
        if x >= 0:
            sign = "+"
        else:
            sign = "-"
        
        # Output the full equation for each energy level
        print(f"E{i + 1} = α {sign} {abs(round(x, 4))}β")

# Run the calculation and print the results
calculate_glyoxal_energies()