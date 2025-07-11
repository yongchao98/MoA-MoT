import numpy as np

def solve_glyoxal_huckel():
    """
    Calculates the 4 pi-electron energy levels of glyoxal using modified Huckel theory.
    """
    # Parameters for Oxygen from the problem description
    h_O = 1.0
    k_CO = 0.8
    
    # Standard parameters for Carbon
    h_C = 0.0
    k_CC = 1.0

    # Construct the Huckel matrix for glyoxal O(1)-C(2)-C(3)-O(4)
    # The elements are derived from h_x and k_xy as explained in the plan.
    huckel_matrix = np.array([
        [h_O,  k_CO, 0.0,  0.0 ],
        [k_CO, h_C,  k_CC, 0.0 ],
        [0.0,  k_CC, h_C,  k_CO],
        [0.0,  0.0,  k_CO, h_O ]
    ])

    # Calculate the eigenvalues of the Huckel matrix
    # These eigenvalues are the 'x' values in the energy equation E = α + xβ
    eigenvalues = np.linalg.eigvals(huckel_matrix)

    # The resonance integral β is negative, so a larger eigenvalue 'x' corresponds
    # to a lower (more stable) energy level. We sort them in descending order.
    sorted_eigenvalues = np.sort(eigenvalues)[::-1]

    # Print the final energy equations
    print("The 4 pi-electron energies (E) for glyoxal are found from the eigenvalues (x) of the Huckel matrix,")
    print("where E = α + xβ.")
    print("Since the resonance integral β is negative, the lowest energy corresponds to the largest x.\n")

    print("The energies are:")
    for i, x in enumerate(sorted_eigenvalues):
        if x >= 0:
            sign = "+"
            value = x
        else:
            sign = "-"
            value = -x  # Use the absolute value for printing

        # The problem asks to output the equation for each energy level
        print(f"E{i+1} = α {sign} {value:.4f}β")

# Execute the function to find and print the energies
solve_glyoxal_huckel()