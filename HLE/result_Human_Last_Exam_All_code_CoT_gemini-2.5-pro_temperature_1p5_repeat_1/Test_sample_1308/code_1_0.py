import numpy as np

def calculate_glyoxal_energies():
    """
    Calculates the 4 pi-electron molecular orbital energies of glyoxal
    using modified Huckel theory.
    """
    # Parameters for the heteroatom (Oxygen)
    h_O = 1.0
    k_CO = 0.8

    # The Huckel matrix is set up for the O1-C2-C3-O4 system.
    # The elements are coefficients for beta, derived from H' = (H - αI) / β.
    # Diagonal elements are h_x values (h_C=0, h_O=1).
    # Off-diagonal elements are k_xy values (k_CC=1, k_CO=0.8).
    H_prime = np.array([
        [h_O, k_CO,  0,    0],
        [k_CO,  0,    1,    0],
        [0,     1,    0,  k_CO],
        [0,     0,  k_CO, h_O]
    ])

    # The eigenvalues of this matrix correspond to the 'x' values in E = α + xβ
    eigenvalues = np.linalg.eigvals(H_prime)

    # Since β is negative, a larger eigenvalue corresponds to a lower (more stable) energy.
    # We sort the eigenvalues in descending order to list energies from lowest to highest.
    sorted_eigenvalues = np.sort(eigenvalues)[::-1]

    print("The 4 Huckel molecular orbital energies for glyoxal are:")
    # We will now print each of the 4 energies in the final equation form
    for i, x in enumerate(sorted_eigenvalues):
        # Format the output to show the equation clearly.
        # This includes showing each number (the eigenvalue) in the equation.
        sign = '+' if x >= 0 else '-'
        print(f"E_{i+1} = α {sign} {abs(x):.4f}β")
        
    # As an additional calculation, for a 4 pi-electron system like glyoxal,
    # the total pi energy is calculated by filling the two lowest energy orbitals (E1, E2)
    # with two electrons each.
    total_pi_energy_coeff = 2 * sorted_eigenvalues[0] + 2 * sorted_eigenvalues[1]
    print("\nFor a 4-electron system, the total π-electron energy is:")
    print(f"E_π = 2*E_1 + 2*E_2 = 4α + {total_pi_energy_coeff:.4f}β")


calculate_glyoxal_energies()