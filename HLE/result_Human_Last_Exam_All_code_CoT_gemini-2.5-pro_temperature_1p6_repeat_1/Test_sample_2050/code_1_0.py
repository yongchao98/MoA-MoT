import math

def calculate_ground_state_energy():
    """
    Calculates and prints the ground state energy for the given system.
    """
    # Parameters from the problem description
    N = 7  # Number of cavities
    M = 4  # Number of photons

    print("This script calculates the ground state energy E_g for a system of N coupled cavities and M photons.")
    print("The Hamiltonian is H = \u03A3_i [\u03C9*a_i\u2020a_i + U/2*a_i\u2020a_i\u2020a_i a_i - J*(a_i\u2020a_{i+1} + a_i\u2020a_{i-1})].")
    print(f"In the limit U -> \u221E with N={N} and M={M} photons, the ground state energy can be calculated.")
    
    print("\nThe analytical formula for the ground state energy is derived by summing the M lowest single-particle energies of the hopping Hamiltonian:")
    
    # The M=4 lowest single-particle modes correspond to quantum numbers k = 0, 1, 6, and 2.
    # The corresponding sum of cosines for the hopping energy term is:
    # cos(2*pi*0/7) + cos(2*pi*1/7) + cos(2*pi*6/7) + cos(2*pi*2/7)
    # Since cos(2*pi*6/7) = cos(2*pi/7), this simplifies to 1 + 2*cos(2*pi/7) + cos(4*pi/7)
    
    # Each term is a number in the final equation.
    term_omega = M
    term_hopping_prefactor = -2
    sum_cos_term_1 = 1 # for k=0
    sum_cos_term_2_coeff = 2 # for k=1 and k=6
    sum_cos_term_2_angle_num = 2 # for 2*pi
    sum_cos_term_2_angle_den = N # for N
    sum_cos_term_3_coeff = 1 # for k=2
    sum_cos_term_3_angle_num = 4 # for 4*pi
    sum_cos_term_3_angle_den = N # for N
    
    print(f"E_g = {term_omega}*\u03C9 + ({term_hopping_prefactor})*J * ( {sum_cos_term_1} + {sum_cos_term_2_coeff}*cos({sum_cos_term_2_angle_num}*\u03C0/{sum_cos_term_2_angle_den}) + {sum_cos_term_3_coeff}*cos({sum_cos_term_3_angle_num}*\u03C0/{sum_cos_term_3_angle_den}) )")
    
    # Numerical calculation
    cos1 = math.cos(sum_cos_term_2_angle_num * math.pi / sum_cos_term_2_angle_den)
    cos2 = math.cos(sum_cos_term_3_angle_num * math.pi / sum_cos_term_3_angle_den)
    hopping_coeff = term_hopping_prefactor * (sum_cos_term_1 + sum_cos_term_2_coeff * cos1 + sum_cos_term_3_coeff * cos2)

    print("\nEvaluating the expression numerically gives:")
    print(f"E_g = {M}*\u03C9 + ({hopping_coeff:.8f})*J")

calculate_ground_state_energy()