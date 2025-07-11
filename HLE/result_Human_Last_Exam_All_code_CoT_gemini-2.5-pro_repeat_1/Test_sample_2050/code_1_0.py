import numpy as np

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for a system of hard-core bosons on a ring.

    The problem involves N=7 cavities and P=4 photons in the U -> infinity limit.
    This corresponds to the ground state of 4 non-interacting spinless fermions
    on a 7-site tight-binding ring.
    The Hamiltonian is H = -J * sum(a_i^dag a_{i+1} + h.c.).
    The term proportional to omega is a constant offset and is ignored.
    """

    # System parameters
    N = 7  # Number of cavities (sites)
    P = 4  # Number of photons (particles)
    
    # We calculate the energy in units of J
    J = 1.0

    # The single-particle energy eigenvalues of a tight-binding ring are E_m = -2J*cos(2*pi*m/N)
    m_values = np.arange(N)
    k_values = 2 * np.pi * m_values / N
    single_particle_energies = -2 * J * np.cos(k_values)

    # The ground state energy of the P-particle system is the sum of the P lowest single-particle energies.
    sorted_energies = np.sort(single_particle_energies)
    ground_state_energy = np.sum(sorted_energies[:P])

    # --- Output the explanation and results ---
    print("The problem considers a ring of N=7 cavities with P=4 photons, in the limit of infinite on-site repulsion U.")
    print("This system maps to 4 non-interacting hard-core bosons (or fermions in 1D) on a 7-site ring.")
    print("The ground state energy is the sum of the lowest 4 single-particle energy levels of the tight-binding model.\n")
    print(f"The single-particle energy levels E_m = -2*J*cos(2*pi*m/{N}) for m = 0..{N-1} are:")
    print([f"{e:.6f}*J" for e in sorted_energies])
    
    # The occupied states for the ground state correspond to m=0, m=1, m=6 (degenerate), and one of m=2, m=5 (degenerate).
    # E_gs = E_0 + E_1 + E_6 + E_2
    # E_gs = -2*J*(cos(0) + cos(2pi/7) + cos(12pi/7) + cos(4pi/7))
    # E_gs = -2*J*(1 + 2*cos(2pi/7) + cos(4pi/7))
    
    c1 = np.cos(2*np.pi/N)
    c2 = np.cos(4*np.pi/N)
    sum_of_cos = 1 + 2*c1 + c2

    print("\nThe ground state energy is the sum of the four lowest levels:")
    print(f"E_gs = ( {sorted_energies[0]:.6f} + {sorted_energies[1]:.6f} + {sorted_energies[2]:.6f} + {sorted_energies[3]:.6f} ) * J")
    print(f"E_gs = {ground_state_energy:.6f} * J\n")
    
    print("This can also be expressed analytically:")
    print("E_gs = -2 * J * (1 + 2*cos(2*pi/7) + cos(4*pi/7))")
    print("Substituting the values:")
    print(f"E_gs = -2 * J * (1 + 2*({c1:.6f}) + ({c2:.6f}))")
    print(f"E_gs = -2 * J * ({sum_of_cos:.6f})")
    print(f"E_gs = {ground_state_energy:.6f} * J")
    print("\nNote: The total energy has an additional term 4*omega, which is a constant offset and has been set to 0.")

calculate_ground_state_energy()