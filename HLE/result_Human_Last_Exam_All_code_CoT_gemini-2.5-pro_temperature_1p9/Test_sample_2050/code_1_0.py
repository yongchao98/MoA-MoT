import math

def solve_ground_state_energy():
    """
    Calculates the ground state energy for 4 hard-core bosons on a 7-site ring.
    """
    N = 7  # Number of cavities
    P = 4  # Number of photons

    # The first term in the Hamiltonian is omega * n_total.
    # Since the number of photons is fixed at P, this contributes P * omega to the energy.
    c_omega = P

    # The second term is the hopping term. In the hard-core limit, we find the ground state
    # by filling the lowest single-particle energy levels of a tight-binding ring.
    # The single-particle energies are E_m = -2 * J * cos(2 * pi * m / N).
    # We need to find the P=4 lowest energy states.

    # Possible m values for N=7, and their degeneracies:
    # m=0 (1 state)
    # m=+/-1 (2 states)
    # m=+/-2 (2 states)
    # m=+/-3 (2 states)
    
    # We rank the energies: E_0 < E_{ +/-1} < E_{ +/-2} < E_{ +/-3}
    # To find the ground state for 4 particles, we fill the 4 lowest energy levels:
    # 1. m=0 state
    # 2. m=1 state
    # 3. m=-1 state
    # 4. one of the m=+/-2 states (they are degenerate, so we take m=2)
    
    occupied_m_values = [0, 1, -1, 2]
    
    # The coefficient of J is the sum of -2 * cos(2 * pi * m / N) for these m values.
    c_j = 0
    print("The ground state energy E_ground is the sum of a constant energy term and the hopping energies.")
    print(f"The constant energy is {c_omega} * omega.")
    print("The hopping energy is E_J = J * Sum[-2 * cos(2*pi*m/N)] for the lowest 4 energy states.")
    print("\nThe occupied single-particle states are for m =", occupied_m_values)
    
    print("\nThe energy contribution for each occupied state (in units of J):")
    for m in occupied_m_values:
        energy_contribution = -2 * math.cos(2 * math.pi * m / N)
        print(f"  m = {m:>2}: {-2 * math.cos(2 * math.pi * m / N):.8f} * J")
        c_j += energy_contribution

    # Print the final equation for the ground state energy
    print("\nThe final ground state energy is E_ground = C_omega * omega + C_J * J")
    print(f"The coefficient for omega is the number of photons:")
    print(f"C_omega = {c_omega}")
    print("\nThe coefficient for J is the sum of the individual state energies:")
    print(f"C_J = {c_j:.8f}")

solve_ground_state_energy()