import numpy as np

def solve_statistical_mechanics():
    """
    Solves the statistical mechanics problem to find the partition function
    and the most probable number of particles in each energy level.
    """
    # Step 1: Define Constants
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T = 4200  # Kelvin
    k_B = 8.617333262e-5  # Boltzmann constant in eV/K
    total_moles = 3.0 # N = 3 * N_A

    # The thermal energy scale
    kT = k_B * T

    n_values = np.arange(1, 6)
    
    # Step 2: Calculate Energy Levels (E_n)
    # Calculate epsilon_n for n=1 to 5
    epsilon_n = epsilon_0_eV * (1 + n_values**n_values)
    
    # Calculate cumulative energy levels E_n
    E_n = np.cumsum(epsilon_n)

    # Step 3: Calculate Degeneracies (g_n)
    g_n = 2 * n_values + 1

    # Step 4: Calculate the Single-Particle Partition Function (Z)
    # Calculate the Boltzmann factor for each level
    boltzmann_factors = np.exp(-E_n / kT)
    
    # Calculate each term in the partition function sum
    Z_terms = g_n * boltzmann_factors
    
    # The single-particle partition function is the sum of these terms
    Z = np.sum(Z_terms)
    
    print(f"The single-particle partition function Z is: {Z:.4f}\n")

    # Step 5 & 6: Calculate the most probable number of moles in each energy level
    # The probability of a particle being in level n is P_n = (term_n) / Z
    # The number of moles in level n is moles_n = total_moles * P_n
    moles_in_levels = total_moles * Z_terms / Z
    
    # Step 7: Format the Output
    print("The most probable number of particles in each energy level (in moles) is:")
    result_tuple = tuple(moles_in_levels)
    print(f"({result_tuple[0]:.4f}, {result_tuple[1]:.4f}, {result_tuple[2]:.4f}, {result_tuple[3]:.4f}, {result_tuple[4]:.4f})")
    
    # Final answer in the required format
    final_answer = f"({result_tuple[0]:.4f}, {result_tuple[1]:.4f}, {result_tuple[2]:.4f}, {result_tuple[3]:.4f}, {result_tuple[4]:.4f})"
    print(f"\n<<<{final_answer}>>>")

solve_statistical_mechanics()