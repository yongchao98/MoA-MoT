import numpy as np

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for a system of hard-core bosons on a ring.

    The Hamiltonian is given by:
    H = Sum_i (w * a_i^dag * a_i + 0.5 * U * a_i^dag * a_i^dag * a_i * a_i - J * a_i^dag * a_{i+1} - J * a_i^dag * a_{i-1})
    with N=7 cavities and 4 photons.
    In the limit U -> infinity, no two photons can occupy the same site.
    This corresponds to a system of non-interacting fermions on a ring.

    The single-particle energy eigenvalues of the hopping Hamiltonian (-J * Sum_i (c_i^dag*c_{i+1} + h.c.))
    on a ring of N sites are E_m = -2 * J * cos(2 * pi * m / N).

    The total ground state energy is the sum of the energies of the P lowest occupied levels.
    """
    # System parameters
    N = 7  # Number of cavities
    P = 4  # Number of photons

    # Generate all N single-particle energy eigenvalues (in units of J)
    # E_m = -2 * cos(2 * pi * m / N)
    energies = []
    for m in range(N):
        energies.append(-2 * np.cos(2 * np.pi * m / N))

    # Sort the energies in ascending order
    energies.sort()

    # The ground state energy is the sum of the lowest P energies
    lowest_p_energies = energies[0:P]
    ground_state_energy_J = sum(lowest_p_energies)

    # Print the explanation and the result as an equation
    print(f"For N={N} sites and P={P} photons in the hard-core limit (U -> infinity),")
    print("we fill the P lowest single-particle energy states of the tight-binding ring.")
    
    print("\nThe single-particle energies sorted from lowest to highest are (in units of J):")
    # Print all energies for context
    for i, e in enumerate(energies):
        print(f"Level {i+1}: {e:.4f} J")

    print(f"\nThe ground state energy is the sum of the {P} lowest energies:")
    
    energy_terms = [f"{e:.4f}" for e in lowest_p_energies]
    # Build the string for the equation
    # Adding parentheses for each number to make it clear, especially for negative numbers
    equation_terms_str = " + ".join([f"({term})" for term in energy_terms])
    
    print(f"E_ground = [ {equation_terms_str} ] * J")
    print(f"E_ground = {ground_state_energy_J:.4f} * J")
    print("\nNote: A constant energy offset of 4*omega from the on-site energy term is not included in this result.")
    
    return ground_state_energy_J

if __name__ == '__main__':
    final_energy_coeff = calculate_ground_state_energy()
    # The final answer format is specific, returning the numerical coefficient of J
    # print(f"<<<{final_energy_coeff:.4f}>>>")

calculate_ground_state_energy()