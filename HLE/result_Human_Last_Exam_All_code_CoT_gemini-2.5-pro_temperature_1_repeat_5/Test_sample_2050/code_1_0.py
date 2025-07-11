import numpy as np

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for a system of hard-core bosons on a ring.
    """
    # System parameters from the problem description
    N = 7  # Number of cavities
    M = 4  # Number of photons

    # In the U -> infinity limit, we have hard-core bosons.
    # The problem reduces to finding the ground state of M non-interacting
    # particles on an N-site ring. The Hamiltonian for the hopping part is
    # H_hop = -J * sum(a_i^+ a_{i+1} + h.c.).
    # The single-particle energy levels of H_hop are E_m = -2*J*cos(2*pi*m/N).

    # We will calculate the energy coefficients (E_m / J) for each momentum state m.
    energy_levels = []
    for m in range(N):
        # The energy coefficient is -2 * cos(2*pi*m/N)
        coeff = -2 * np.cos(2 * np.pi * m / N)
        energy_levels.append((coeff, m))

    # Sort the energy levels in ascending order. The first element of the tuple (the energy)
    # is used for sorting.
    energy_levels.sort()

    # The ground state is formed by filling the M lowest energy levels.
    occupied_levels = energy_levels[:M]

    # The total hopping energy is the sum of the energies of these M levels.
    total_hopping_coeff = sum(level[0] for level in occupied_levels)
    
    # --- Output the results step-by-step ---
    
    print(f"The system consists of N={N} cavities and M={M} photons in a ring.")
    print(f"In the limit U -> infinity, the photons act as hard-core bosons, and the problem simplifies.")
    print(f"The ground state energy is the sum of the {M} lowest single-particle energy levels from the hopping Hamiltonian.")
    
    occupied_m_values = [m for _, m in occupied_levels]
    print(f"\nThe {M} lowest single-particle energy levels correspond to quantum numbers m = {occupied_m_values}.")

    print("\nThe contributions to the hopping energy from these levels are:")
    
    equation_parts = []
    for energy_coeff, m in occupied_levels:
        print(f"  For m={m}: E_{m}/J = {energy_coeff:.4f}")
        equation_parts.append(f"({energy_coeff:.4f})")

    # The total ground state energy is E_ground = M*omega + E_hopping
    print("\nThe total ground state energy E_ground is given by the expression:")
    
    # We print the equation with each number explicitly, as requested.
    final_equation_str = f"E_ground = {M}*omega + [" + " + ".join(equation_parts) + "] * J"
    print(final_equation_str)
    
    # We also print the final simplified expression.
    final_simplified_str = f"E_ground = {M}*omega + {total_hopping_coeff:.4f}*J"
    print("\nWhich simplifies to:")
    print(final_simplified_str)


if __name__ == '__main__':
    calculate_ground_state_energy()