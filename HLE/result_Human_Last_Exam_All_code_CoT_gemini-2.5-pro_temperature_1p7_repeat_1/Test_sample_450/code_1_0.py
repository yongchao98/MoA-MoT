import math
import numpy as np

def solve_particle_distribution():
    """
    Solves for the most probable number of particles in each energy level for a given system.
    """
    # Step 1: Define Constants
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T = 4200.0  # Kelvin
    k_B_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K
    total_moles = 3.0

    # Number of energy levels
    num_levels = 5
    n_values = range(1, num_levels + 1)

    # Step 2: Calculate Energy Levels and Degeneracies
    g_n = [2 * n + 1 for n in n_values]
    
    # Calculate intermediate epsilon values
    epsilon_i_values = [epsilon_0_eV * (1 + float(n)**n) for n in n_values]
    
    # Calculate cumulative energy levels E_n
    E_n = np.cumsum(epsilon_i_values).tolist()

    # Step 3: Calculate the Single-Particle Partition Function (Z1)
    # The total partition function for N distinguishable particles is Z = (Z1)^N.
    # Calculating Z is computationally infeasible and unnecessary for finding populations.
    # We only need Z1, the single-particle partition function.
    
    beta = 1.0 / (k_B_eV_K * T)
    
    boltzmann_terms = []
    Z1 = 0.0
    for i in range(num_levels):
        term = g_n[i] * math.exp(-E_n[i] * beta)
        boltzmann_terms.append(term)
        Z1 += term
        
    print(f"The single-particle partition function Z1 is: {Z1:.5f}")

    # Step 4: Calculate the Most Probable Number of Moles in Each Level
    moles_in_levels = []
    for i in range(num_levels):
        # Probability P_n = term_n / Z1
        # Number of moles moles_n = total_moles * P_n
        moles = total_moles * boltzmann_terms[i] / Z1
        moles_in_levels.append(moles)
    
    # Step 5: Output the Final Result
    # "output each number in the final equation!" is interpreted as printing the final tuple.
    print("The most probable number of moles in each energy level (E1, E2, E3, E4, E5) is:")
    
    # We construct the final string tuple to be printed
    final_answer_tuple_str = "(" + ", ".join([f"{mole:.5f}" for mole in moles_in_levels]) + ")"
    print(final_answer_tuple_str)
    
    # Present the final answer in the requested format
    final_answer_tuple_obj = tuple(moles_in_levels)
    print(f"\n<<< {final_answer_tuple_obj} >>>")

solve_particle_distribution()