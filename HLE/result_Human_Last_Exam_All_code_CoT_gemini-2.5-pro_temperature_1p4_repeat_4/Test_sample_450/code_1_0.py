import numpy as np

def solve_statistical_mechanics():
    """
    Calculates the partition function and the number of moles in each energy level
    for a system of distinguishable particles.
    """
    # 1. Define Constants
    N_moles = 3.0  # Total number of moles
    epsilon_0 = 6.9e-3  # Base energy in eV
    T = 4200.0  # Temperature in Kelvin
    k_B = 8.617333262e-5  # Boltzmann constant in eV/K

    # Thermal energy
    k_B_T = k_B * T

    # Define the quantum number range
    n = np.arange(1, 6)

    # 2. Calculate Energy Levels
    # Calculate intermediate epsilon_n values
    epsilon_n = epsilon_0 * (1 + n**n)
    # Calculate the actual energy levels E_n by cumulative sum
    E_n = np.cumsum(epsilon_n)

    # 3. Calculate Degeneracies
    g_n = 2 * n + 1

    # 4. Calculate the Single-Particle Partition Function (Z1)
    # Calculate Boltzmann factor terms
    boltzmann_terms = g_n * np.exp(-E_n / k_B_T)
    # Sum the terms to get Z1
    Z1 = np.sum(boltzmann_terms)

    print(f"The single-particle partition function Z1 is: {Z1}")

    # 5. Calculate Occupation Numbers in Moles
    # Calculate the probability of a particle being in each state
    P_n = boltzmann_terms / Z1
    # Calculate the number of moles in each energy level
    moles_n = N_moles * P_n

    print("\nThe most probable number of particles in each energy level (in moles) is:")
    # Print the final ordered set
    # The format "(E1, E2, E3, E4, E5)" in the prompt refers to the values associated
    # with these energy levels.
    print(f"({moles_n[0]}, {moles_n[1]}, {moles_n[2]}, {moles_n[3]}, {moles_n[4]})")
    
    # Format for the final answer block
    final_answer = tuple(moles_n)
    return final_answer

# Execute the function and print the final answer in the required format
result = solve_statistical_mechanics()
print(f'<<<{result}>>>')
