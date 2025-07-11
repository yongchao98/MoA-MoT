import math

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for a system of hard-core bosons on a ring.

    The problem involves N=7 cavities in a ring with M=4 photons, in the limit
    of infinite on-site repulsion (U -> infinity). This corresponds to a system
    of hard-core bosons, which is equivalent to non-interacting spinless fermions.
    The Hamiltonian simplifies to a tight-binding model.
    """
    # Define problem parameters
    N = 7  # Number of cavities
    M = 4  # Number of photons

    # The on-site energy term, Sum(omega * a_i^dagger * a_i), contributes M*omega
    # to the total energy, as the total number of photons is fixed at M.
    omega_coefficient = M

    # The hopping part of the Hamiltonian, -J * Sum(a_i^dagger * a_{i+1} + h.c.),
    # is solved by diagonalizing it in momentum space. The single-particle
    # energy eigenstates have energies E_m = -2*J*cos(k_m), where k_m = 2*pi*m/N.

    # To find the ground state, we fill the M=4 lowest energy levels.
    # The lowest energies correspond to the largest values of cos(2*pi*m/N).
    # The allowed integer momentum quantum numbers for a ring of size N=7 are
    # m = 0, +/-1, +/-2, +/-3.
    all_m = list(range(-N // 2, N // 2 + 1))

    # Create a dictionary of m -> cos(2*pi*m/N) to represent the energy ordering term.
    energies = {m: math.cos(2 * math.pi * m / N) for m in all_m}

    # Sort the m values by their corresponding cos value in descending order
    # to find the states with the lowest energies.
    sorted_m = sorted(energies, key=energies.get, reverse=True)

    # The 4 photons will occupy the M lowest-energy states.
    occupied_m_states = sorted_m[:M]

    # The total energy from the hopping term is the sum of the energies of
    # these occupied states: E_hopping = Sum over occupied m of -2*J*cos(2*pi*m/N).
    # We calculate the numeric coefficient for J.
    sum_of_cosines = sum(energies[m] for m in occupied_m_states)
    J_coefficient = -2 * sum_of_cosines

    # The final ground state energy is E_g = M*omega + E_hopping.
    # We will print the full equation with the calculated numeric coefficients.
    print("The ground state energy E_g is given by the equation:")
    # We check the sign of the J coefficient to format the equation nicely.
    if J_coefficient < 0:
        print(f"E_g = {omega_coefficient} * omega - {-J_coefficient:.4f} * J")
    else:
        print(f"E_g = {omega_coefficient} * omega + {J_coefficient:.4f} * J")

if __name__ == '__main__':
    calculate_ground_state_energy()
