import numpy as np

def solve_huckel_energies():
    """
    Calculates the pi orbital energies of glyoxal using the modified Huckel theory.
    """

    # Step 1: Define the Huckel parameters from the problem description.
    # For the heteroatom (Oxygen):
    h_O = 1.0
    k_CO = 0.8
    # For carbon atoms, standard parameters are used:
    h_C = 0.0
    k_CC = 1.0

    # Step 2: Construct the Huckel matrix for the O1-C2-C3-O4 pi system.
    huckel_matrix = np.array([
        [h_O, k_CO, 0.0,  0.0],
        [k_CO, h_C, k_CC, 0.0],
        [0.0,  k_CC, h_C, k_CO],
        [0.0,  0.0, k_CO, h_O]
    ])

    # Step 3: Calculate the eigenvalues of the Huckel matrix.
    # The energies E are given by E = α + λβ, where λ are the eigenvalues.
    eigenvalues = np.linalg.eigvals(huckel_matrix)

    # Step 4: Sort the eigenvalues to represent ordered energy levels (from lowest to highest).
    eigenvalues.sort()

    # Step 5: Print the four resulting energies in terms of α (alpha) and β (beta).
    # Each number is explicitly shown in the final equation as requested.
    print("The four pi system energies (E) for glyoxal are:")
    for i, lam in enumerate(eigenvalues):
        print(f"E{i+1} = \u03B1 + ({lam:.4f})\u03B2")

solve_huckel_energies()