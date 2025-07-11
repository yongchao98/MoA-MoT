import numpy as np

def solve_statistical_mechanics():
    """
    Calculates the partition function and the most probable number of particles
    in each energy level for a given system.
    """
    # Step 1: Define Physical Constants and System Parameters
    N_moles = 3.0  # Total number of moles
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T_K = 4200.0  # Temperature in Kelvin
    k_B_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K

    # Thermal energy
    kBT = k_B_eV_K * T_K

    num_levels = 5
    g_n_list = []
    E_n_list = []

    # Step 2: Calculate Energy Levels (E_n) and Degeneracies (g_n)
    E_cumulative_eV = 0.0
    for n in range(1, num_levels + 1):
        # Degeneracy
        g_n = 2 * n + 1
        g_n_list.append(g_n)

        # Energy increment
        epsilon_n = epsilon_0_eV * (1 + n**n)
        
        # Cumulative energy level
        E_cumulative_eV += epsilon_n
        E_n_list.append(E_cumulative_eV)

    # Step 3: Compute the Single-Particle Partition Function (Z_1)
    Z_single = 0.0
    boltzmann_terms = []
    for i in range(num_levels):
        g_n = g_n_list[i]
        E_n = E_n_list[i]
        
        # Boltzmann factor term for the partition function sum
        term = g_n * np.exp(-E_n / kBT)
        boltzmann_terms.append(term)
        Z_single += term

    print(f"The single-particle partition function Z_1 is: {Z_single:.5f}")

    # Step 4: Calculate the Number of Moles in Each Energy Level
    moles_in_levels = []
    for i in range(num_levels):
        # Probability of a particle being in level n
        prob_n = boltzmann_terms[i] / Z_single
        
        # Number of moles in level n
        moles_n = N_moles * prob_n
        moles_in_levels.append(moles_n)
    
    # Step 5: Format and Print the Final Output
    # The problem asks to output each number in the final equation.
    # We will print the tuple of mole numbers.
    final_moles_tuple = tuple(moles_in_levels)
    print("\nThe most probable number of particles in each energy level (in moles) is:")
    print(f"(moles_E1, moles_E2, moles_E3, moles_E4, moles_E5) = {final_moles_tuple}")
    
    # Final answer in the required format
    answer = f"<<<{final_moles_tuple}>>>"
    # This print is for the platform, but the user sees the one above.
    # print(answer)

solve_statistical_mechanics()