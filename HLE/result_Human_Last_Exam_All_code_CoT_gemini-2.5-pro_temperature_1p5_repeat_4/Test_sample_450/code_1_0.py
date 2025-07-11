import numpy as np
from scipy import constants

def calculate_moles_distribution():
    """
    Calculates the partition function and the most probable number of particles
    in each energy level for a given statistical mechanics system.
    """
    # 1. Define Constants and Parameters
    epsilon_0_meV = 6.9  # meV
    T_kelvin = 4200.0    # Kelvin
    total_moles = 3.0    # Total moles is 3 (from N = 3*N_A)

    # Physical constants
    k_B = constants.k      # Boltzmann constant in J/K
    e_charge = constants.e # Elementary charge in C

    # Convert epsilon_0 from meV to Joules for consistency
    epsilon_0_J = epsilon_0_meV * 1e-3 * e_charge

    # Thermal energy factor for the exponent
    beta = 1.0 / (k_B * T_kelvin)

    # 2. Calculate Energy Levels and Degeneracies for n=1 to 5
    n_values = np.arange(1, 6)

    # Degeneracies: g_n = 2n + 1
    g_n_values = 2 * n_values + 1

    # Energy increments: epsilon_n = epsilon_0 * (1 + n^n)
    # Use float type for n to handle n^n calculation robustly
    epsilon_n_values = epsilon_0_J * (1 + n_values.astype(float)**n_values.astype(float))

    # Cumulative energy levels: E_n = sum of epsilon_i up to n
    E_n_values = np.cumsum(epsilon_n_values)

    # 3. Calculate the Single-Particle Partition Function (Z)
    boltzmann_terms = g_n_values * np.exp(-E_n_values * beta)
    Z = np.sum(boltzmann_terms)

    # 4. Calculate Population in Moles
    # Probability for a particle to be in level n
    probabilities = boltzmann_terms / Z
    
    # Most probable number of moles in level n
    moles_n = total_moles * probabilities

    # 5. Output the Results
    print(f"The single-particle partition function Z is: {Z:.6f}")
    
    print("\nThe most probable number of moles for each energy level (the 'final equation'):")
    total_moles_check = 0
    for i in range(len(n_values)):
        print(f"Number of moles in E_{i+1}: {moles_n[i]:.6f}")
        total_moles_check += moles_n[i]

    print(f"\nSum of moles: {total_moles_check:.6f} (should be close to {total_moles})")
    
    final_ordered_set = tuple(moles_n)
    print(f"\nThe ordered set of moles (E1, E2, E3, E4, E5) is: {final_ordered_set}")

    # Final answer formatted as requested
    print(f"\n<<<{final_ordered_set}>>>")


if __name__ == "__main__":
    calculate_moles_distribution()