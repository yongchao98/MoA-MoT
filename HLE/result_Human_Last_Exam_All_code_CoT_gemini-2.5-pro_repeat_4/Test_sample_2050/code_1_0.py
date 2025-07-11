import math

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for 4 hard-core bosons on a 7-site ring.

    The ground state energy is E_g = 4*omega + E_hopping.
    E_hopping is the sum of the 4 lowest single-particle energies.
    Single-particle energies: epsilon_m = -2*J*cos(2*pi*m/N) with N=7.
    The 4 lowest energy states correspond to m = 0, +1, -1, +2.
    """
    N = 7

    # Calculate the values of the cosine terms
    cos_0 = math.cos(0)
    cos_2pi_over_7 = math.cos(2 * math.pi / N)
    cos_4pi_over_7 = math.cos(4 * math.pi / N)

    # Calculate the energy contribution from each occupied mode in terms of J
    # E_m = -2J * cos(2*pi*m/N)
    e0_coeff = -2 * cos_0
    e1_and_minus1_coeff = 2 * (-2 * cos_2pi_over_7)
    e2_coeff = -2 * cos_4pi_over_7

    # Total coefficient for J in the hopping energy
    total_j_coeff = e0_coeff + e1_and_minus1_coeff + e2_coeff

    # Print the step-by-step calculation
    print("The ground state energy E_g is the sum of the on-site energy and the hopping energy:")
    print("E_g = 4*omega + E_hopping")
    print("\nThe hopping energy is the sum of the energies of the 4 lowest occupied single-particle modes (m=0, +/-1, +2):")
    print(f"E_hopping = (E_0 + E_1 + E_-1 + E_2)")
    print(f"E_hopping = (-2*cos(0) - 4*cos(2*pi/7) - 2*cos(4*pi/7)) * J")
    print("\nLet's calculate the numerical values of the coefficients of J:")
    print(f"Contribution from m=0: -2*cos(0) = {e0_coeff:.5f}")
    print(f"Contribution from m=+/-1: -4*cos(2*pi/7) = {e1_and_minus1_coeff:.5f}")
    print(f"Contribution from m=+2: -2*cos(4*pi/7) = {e2_coeff:.5f}")
    
    print("\nThe final equation for the ground state energy is:")
    print(f"E_g = 4*omega + ({e0_coeff:.5f} + ({e1_and_minus1_coeff:.5f}) + {e2_coeff:.5f}) * J")
    print(f"E_g = 4*omega + {total_j_coeff:.5f} * J")

calculate_ground_state_energy()
