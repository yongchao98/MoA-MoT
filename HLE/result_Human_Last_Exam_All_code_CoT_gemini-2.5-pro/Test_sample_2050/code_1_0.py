import math

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for 4 hard-core bosons on a 7-site ring.
    
    The Hamiltonian is H = Sigma_i(omega*n_i + U/2*n_i*(n_i-1) - J*(a_i^dag a_{i+1} + h.c.)).
    In the limit U -> infinity, n_i can only be 0 or 1.
    The problem maps to 4 non-interacting fermions on a 7-site ring.
    The ground state energy is the sum of the on-site energy (4*omega) and the
    hopping energy, which is the sum of the 4 lowest single-particle energies.
    
    Single-particle energies: E_m = -2*J*cos(2*pi*m/N)
    For N=7, the lowest 4 levels correspond to m=0, m=+1, m=-1, and m=+2 (or -2).
    """
    N = 7
    num_photons = 4

    print("Step 1: Understand the model in the U -> infinity limit.")
    print("The U -> infinity term forces each cavity to have at most one photon.")
    print(f"Thus, with {num_photons} photons and {N} cavities, we must place 4 photons on 4 different sites.")
    print("The problem becomes equivalent to 4 non-interacting fermions on a ring.")
    print("\nStep 2: Find the single-particle energy levels.")
    print("The single-particle energies for a particle on a ring of N sites are:")
    print("E_m = -2 * J * cos(2 * pi * m / N), where m = 0, +/-1, +/-2, ...")
    
    print("\nStep 3: Find the ground state energy by filling the lowest 4 energy levels.")
    print("For N=7, the lowest energy states correspond to m=0, m=+1, m=-1, and m=+2.")
    print("The total hopping energy is the sum of the energies of these states:")
    print("E_hopping = E_0 + E_1 + E_{-1} + E_2")
    print("E_hopping = [-2J*cos(0)] + [-2J*cos(2*pi/7)] + [-2J*cos(-2*pi/7)] + [-2J*cos(4*pi/7)]")
    
    # Calculate the numerical values for the terms in the equation
    # The term for m=0 gives cos(0) = 1
    term_m0 = 1.0
    # The terms for m=+1 and m=-1 are degenerate, cos(x) = cos(-x)
    term_m1 = 2 * math.cos(2 * math.pi / N)
    # The term for m=+2
    term_m2 = math.cos(4 * math.pi / N)
    
    # The total coefficient for -2J
    sum_of_cosines = term_m0 + term_m1 + term_m2
    # The final coefficient for J
    j_coefficient = -2 * sum_of_cosines
    
    print("\nStep 4: Calculate the numerical value of the hopping energy.")
    print("The total ground state energy is E_g = 4*omega + E_hopping.")
    print("Let's calculate E_hopping. The final equation is:")
    print(f"E_hopping = -2*J * (cos(2*pi*0/7) + 2*cos(2*pi*1/7) + cos(2*pi*2/7))")
    print(f"E_hopping = -2*J * ( {term_m0:.5f} + {term_m1:.5f} + {term_m2:.5f} )")
    print(f"E_hopping = -2*J * ( {sum_of_cosines:.5f} )")
    print(f"E_hopping = ( {j_coefficient:.5f} ) * J")

    print("\n-------------------------------------------------------------")
    print("Final Answer: The ground state energy of the system is:")
    print(f"E_g = 4*omega + {j_coefficient:.5f}*J")
    print("-------------------------------------------------------------")

if __name__ == '__main__':
    calculate_ground_state_energy()