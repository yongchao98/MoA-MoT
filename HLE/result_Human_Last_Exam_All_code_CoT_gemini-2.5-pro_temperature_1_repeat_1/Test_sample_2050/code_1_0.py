import numpy as np

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for N hard-core bosons on a ring.

    The U -> infinity limit forces the bosons to be "hard-core", meaning no two
    bosons can occupy the same site. This system is equivalent to non-interacting
    fermions on a lattice. The ground state is found by filling the lowest
    single-particle energy levels.
    """
    # System parameters
    N = 7
    num_photons = 4

    print(f"System: {N} cavities in a ring with {num_photons} photons.")
    print(f"Analysis in the U -> infinity limit (hard-core bosons).\n")

    # The single-particle energy levels for the hopping term (-J*sum(a_i+ a_{i+1} + h.c.))
    # on a ring are E'_m = -2J*cos(2*pi*m/N). We calculate the coefficient of J.
    m_values = np.arange(N)
    energies_J_coeff = -2 * np.cos(2 * np.pi * m_values / N)

    # To find the ground state, we fill the lowest energy levels.
    # We sort the energies and take the sum of the lowest `num_photons` values.
    sorted_energies_J_coeff = np.sort(energies_J_coeff)
    occupied_energies_J_coeff = sorted_energies_J_coeff[:num_photons]
    total_energy_J_coeff = np.sum(occupied_energies_J_coeff)

    # The total energy is E_gs = N_photons*w + (sum of occupied E'_m)
    print("The total ground state energy is the sum of the constant w-term and the J-hopping term.")
    print("E_gs = (Number of photons)*w + (Sum of lowest single-particle hopping energies)*J\n")
    print("The energies of the 4 occupied levels (coefficient of J) are:")
    for i, energy_coeff in enumerate(occupied_energies_J_coeff):
        print(f"  Level {i+1}: {energy_coeff:.4f}")

    # Display the final equation with each number, as requested.
    sum_string = f"({occupied_energies_J_coeff[0]:.4f})"
    for i in range(1, len(occupied_energies_J_coeff)):
        # Add a plus sign for visual clarity, even though the number is negative
        sum_string += f" + ({occupied_energies_J_coeff[i]:.4f})"

    print(f"\nThe final energy equation is:")
    print(f"E_gs = {num_photons}*w + [ {sum_string} ]*J")
    print(f"E_gs = {num_photons}*w + ({total_energy_J_coeff:.4f})*J")


if __name__ == "__main__":
    calculate_ground_state_energy()
    # The non-trivial numerical part of the answer is the coefficient of J.
    N = 7
    num_photons = 4
    m_values = np.arange(N)
    energies_J_coeff = -2 * np.cos(2 * np.pi * m_values / N)
    sorted_energies_J_coeff = np.sort(energies_J_coeff)
    occupied_energies_J_coeff = sorted_energies_J_coeff[:num_photons]
    total_energy_J_coeff = np.sum(occupied_energies_J_coeff)
    # The final answer is the coefficient of J.
    print(f"\n<<<{total_energy_J_coeff:.4f}>>>")