import numpy as np

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for a system of hard-core bosons on a ring.

    The problem considers N=7 cavities in a ring with K=4 photons, in the
    limit of infinite on-site interaction U. This maps to a system of
    non-interacting spinless fermions.
    """
    N = 7  # Number of cavities
    K = 4  # Number of photons

    # The single-particle energy eigenvalues for a tight-binding ring are E_m = -2*J*cos(2*pi*m/N).
    # We calculate the coefficients of J for all N states.
    m_values = np.arange(N)
    energy_coeffs = -2 * np.cos(2 * np.pi * m_values / N)

    # To find the ground state energy for K fermions, we fill the K lowest
    # available single-particle energy levels. We sort the energies to find them.
    sorted_energy_coeffs = np.sort(energy_coeffs)
    ground_state_coeffs = sorted_energy_coeffs[:K]

    # The total kinetic energy is the sum of the energies of the occupied states.
    total_kinetic_energy_coeff = np.sum(ground_state_coeffs)

    # The total ground state energy is E_gs = K*omega + E_kinetic.
    print("The system maps to 4 non-interacting fermions on a 7-site ring.")
    print("The ground state is found by filling the 4 lowest single-particle energy levels.")
    print("The total energy is E = K*omega + Sum(E_m) for the 4 lowest states m.")
    
    # We print the final equation showing each numerical contribution.
    print("\nThe final equation for the ground state energy E is:")
    
    equation_str = f"E = {K}*omega"
    
    # We build the equation string piece by piece.
    # Note that the values are coefficients of J.
    part_eqs = []
    for coeff in ground_state_coeffs:
        if coeff >= 0:
            part_eqs.append(f"+ {coeff:.4f}*J")
        else:
            part_eqs.append(f"- {-coeff:.4f}*J")
            
    # Python's np.sort might put -1.2470 before or after -1.2470, which is fine
    # because they are degenerate. The printout shows one such combination.
    print(f"E = {K}*omega {part_eqs[0]} {part_eqs[1]} {part_eqs[2]} {part_eqs[3]}")
    
    # Final simplified expression
    print("\nSumming the kinetic terms gives:")
    if total_kinetic_energy_coeff >= 0:
        print(f"E = {K}*omega + {total_kinetic_energy_coeff:.8f}*J")
    else:
        print(f"E = {K}*omega - {-total_kinetic_energy_coeff:.8f}*J")

if __name__ == '__main__':
    calculate_ground_state_energy()