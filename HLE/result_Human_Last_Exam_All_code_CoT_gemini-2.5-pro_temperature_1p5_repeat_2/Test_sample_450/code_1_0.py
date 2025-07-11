import math

def solve_statistical_mechanics():
    """
    Calculates the partition function and the most probable number of moles
    in each energy level for a system of distinguishable particles.
    """
    # 1. Define physical constants and problem parameters
    epsilon_0_eV = 6.9e-3  # Base energy in eV
    T_K = 4200.0           # Temperature in Kelvin
    k_B_eV_K = 8.617333e-5 # Boltzmann constant in eV/K
    N_moles = 3.0          # Total moles of particles
    num_levels = 5

    # 2. Calculate energy levels (E_n) and degeneracies (g_n)
    energies_eV = []
    degeneracies = []
    
    cumulative_epsilon = 0.0
    for n in range(1, num_levels + 1):
        # Degeneracy
        g_n = 2 * n + 1
        degeneracies.append(g_n)
        
        # Energy contribution for level n
        epsilon_n = epsilon_0_eV * (1 + n**n)
        
        # Cumulative energy for the level
        cumulative_epsilon += epsilon_n
        energies_eV.append(cumulative_epsilon)

    # 3. Calculate the single-particle partition function (Z)
    kBT_eV = k_B_eV_K * T_K
    
    Z_terms = []
    for i in range(num_levels):
        g = degeneracies[i]
        E = energies_eV[i]
        # Boltzmann factor term for the level
        term = g * math.exp(-E / kBT_eV)
        Z_terms.append(term)
    
    # The partition function Z is the sum of all terms
    Z = sum(Z_terms)

    # 4. Calculate the most probable number of moles per level
    moles_per_level = []
    for term in Z_terms:
        # Probability P_n = (Boltzmann factor term) / Z
        probability_n = term / Z
        # Moles in level n = Total moles * P_n
        moles_n = N_moles * probability_n
        moles_per_level.append(moles_n)
        
    # 5. Print the results
    print(f"The single-particle partition function Z is: {Z:.5f}")
    
    print("\nThe most probable number of moles in each energy level is given by the ordered set:")
    # The problem requests the final numbers in the "equation"
    # (moles_E1, moles_E2, moles_E3, moles_E4, moles_E5)
    formatted_moles = tuple(moles_per_level)
    print(f"(moles_E1, moles_E2, moles_E3, moles_E4, moles_E5) = {formatted_moles}")


# Execute the function to solve the problem
solve_statistical_mechanics()

# Final answer in the specified format
final_answer = (0.7955938722830883, 1.205434190226305, 0.9895603099905251, 0.009477013898687707, 6.035619175402454e-28)
# We will present this tuple as the final answer
# Let's format it for readability if needed, but the raw tuple is precise.
# For the final <<<...>>> format, a clean tuple representation is best.
final_answer_tuple = (final_answer[0], final_answer[1], final_answer[2], final_answer[3], final_answer[4])
# This line is just for generating the final output and won't be in the block.
# print(f'<<<{final_answer_tuple}>>>')