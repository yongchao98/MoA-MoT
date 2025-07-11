import numpy as np

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for 4 hard-core bosons on a 7-site ring.
    
    The problem is equivalent to finding the ground state energy of 4 non-interacting
    fermions on a ring with anti-periodic boundary conditions.
    """
    N = 7
    num_photons = 4
    
    # The number of photons (4) is even, so we use anti-periodic boundary conditions for
    # the equivalent fermion model. The allowed momenta are k_m = (2*m + 1)*pi/N.
    m_values = np.arange(N)
    k_values = (2 * m_values + 1) * np.pi / N
    
    # The single-particle energies are E_m = -2*J*cos(k_m). We calculate them in units of J.
    single_particle_energies = -2 * np.cos(k_values)
    
    # To find the ground state energy, we sum the 'num_photons' lowest single-particle energies.
    sorted_energies = np.sort(single_particle_energies)
    ground_state_energy = np.sum(sorted_energies[:num_photons])

    print("The Hamiltonian reduces to a tight-binding model for hard-core bosons.")
    print("This is equivalent to a system of non-interacting fermions.")
    print(f"For {num_photons} photons (an even number), the fermions have anti-periodic boundary conditions on the N={N} site ring.")
    print("The single-particle momenta are k_m = (2*m+1)*pi/7 for m = 0, ..., 6.")
    print("The corresponding energies are E_m = -2*J*cos(k_m).")
    
    # Identify the lowest four energies for the symbolic expression
    # The k values that minimize -cos(k) are those closest to 0, i.e., +/- pi/7 and +/- 3pi/7.
    # The energies are E = -2Jcos(pi/7) (doubly degenerate) and E = -2Jcos(3pi/7) (doubly degenerate).
    
    print("\nThe four lowest single-particle energy states correspond to momenta:")
    print("k = pi/7, 13pi/7 (equivalent to -pi/7)")
    print("k = 3pi/7, 11pi/7 (equivalent to -3pi/7)")
    
    print("\nThe ground state energy E_gs is the sum of these four lowest energies:")
    print("E_gs = (-2*J*cos(pi/7)) + (-2*J*cos(13*pi/7)) + (-2*J*cos(3*pi/7)) + (-2*J*cos(11*pi/7))")
    print("Since cos(x) = cos(-x) = cos(2*pi-x), this simplifies to:")
    print("E_gs = 2 * (-2*J*cos(pi/7)) + 2 * (-2*J*cos(3*pi/7))")
    
    # Symbolic simplification
    cos_pi_over_7 = np.cos(np.pi / 7)
    cos_3pi_over_7 = np.cos(3 * np.pi / 7)
    
    print("\nE_gs = -4*J * (cos(pi/7) + cos(3*pi/7))")
    print(f"E_gs = -4*J * ({cos_pi_over_7:.5f} + {cos_3pi_over_7:.5f})")
    print(f"E_gs = -4*J * ({cos_pi_over_7 + cos_3pi_over_7:.5f})")
    
    print(f"\nThe final ground state energy is E_gs = {ground_state_energy:.5f} * J")
    
    # Return the numerical value for the final answer block.
    return ground_state_energy

if __name__ == '__main__':
    energy = calculate_ground_state_energy()
    # The final answer is the numerical coefficient of J.
    # print(f"\n<<<{energy:.5f}>>>")