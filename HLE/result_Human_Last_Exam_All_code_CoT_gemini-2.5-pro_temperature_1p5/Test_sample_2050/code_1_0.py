import numpy as np

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for 4 hard-core bosons on a 7-site ring.
    """
    N = 7  # Number of cavities
    P = 4  # Number of photons (particles)

    # In the U -> infinity limit, the problem maps to non-interacting fermions.
    # The single-particle energy levels for a particle on an N-site ring are
    # given by E_m = -2J * cos(2*pi*m/N). We calculate the coefficient of J.
    
    # We can generate all unique single-particle energies by letting m run from 0 to N-1.
    energy_coeffs = []
    for m in range(N):
        coeff = -2 * np.cos(2 * np.pi * m / N)
        energy_coeffs.append(coeff)

    # To find the ground state, we fill the lowest available energy levels.
    # First, we sort the energies in ascending order.
    energy_coeffs.sort()

    # The ground state for P particles is the sum of the P lowest single-particle energies.
    lowest_four_coeffs = energy_coeffs[:P]
    total_J_coeff = sum(lowest_four_coeffs)

    # The total energy is E = 4*omega + E_hopping.
    # E_hopping is the sum of the energies of the occupied states.
    print("The total ground state energy is E_ground = 4*omega + E_hopping.")
    print("The hopping energy E_hopping is the sum of the four lowest single-particle energy levels.")
    
    print("\nThe contribution of each of the four particles to E_hopping is:")
    term1, term2, term3, term4 = lowest_four_coeffs
    # We output each number in the final equation as requested.
    final_equation_hopping_part = f"E_hopping = ({term1:.4f} + {term2:.4f} + {term3:.4f} + {term4:.4f}) * J"
    print(final_equation_hopping_part)
    
    print(f"\nSumming these values gives the total hopping energy:")
    print(f"E_hopping = {total_J_coeff:.4f} * J")

    print("\nTherefore, the final equation for the ground state energy is:")
    final_equation = f"E_ground = 4*omega + {total_J_coeff:.4f} * J"
    print(final_equation)

calculate_ground_state_energy()