import numpy as np

def solve_statistical_mechanics():
    """
    Calculates the single-particle partition function and the most probable
    number of moles of particles in each energy level for a given system.
    """
    # Step 1: Define constants
    epsilon_0 = 6.9e-3  # eV
    T = 4200.0         # Kelvin
    k_B = 8.617333e-5  # eV/K
    N_total_moles = 3.0 # Total number of particles is 3 * N_A

    num_levels = 5
    
    # Calculate thermal energy
    kbT = k_B * T

    # Step 2: Calculate energy levels (E_n) and degeneracies (g_n)
    n_values = np.arange(1, num_levels + 1)
    
    # Degeneracies g_n = 2n + 1
    g_n = 2 * n_values + 1
    
    # Energy contributions epsilon_n = epsilon_0 * (1 + n^n)
    epsilon_contributions = epsilon_0 * (1 + n_values**n_values)
    
    # Cumulative energy levels E_n = sum_{i=1 to n} epsilon_i
    E_n = np.cumsum(epsilon_contributions)
    
    # Step 3: Calculate single-particle partition function (Z_1)
    # Z_1 = sum(g_n * exp(-E_n / kbT))
    boltzmann_factors = g_n * np.exp(-E_n / kbT)
    Z_1 = np.sum(boltzmann_factors)

    # Step 4: Calculate the most probable number of moles in each level
    # Probability P_n = (g_n * exp(-E_n / kbT)) / Z_1
    # Moles_n = N_total_moles * P_n
    moles_in_level = N_total_moles * boltzmann_factors / Z_1
    
    # Step 5: Output the results
    print(f"The single-particle partition function Z_1 is: {Z_1}")
    print("The most probable number of particles in each energy level (in moles) is:")
    
    # Create the final tuple for the answer
    final_answer = tuple(moles_in_level)
    
    # The final print statement matches the requested format (an ordered set)
    print(f"(E1, E2, E3, E4, E5) = {final_answer}")
    
    # Format the final answer string as requested
    answer_string = f"<<<{final_answer}>>>"
    # This print statement is for direct parsing of the final result.
    # print(answer_string)
    
if __name__ == "__main__":
    solve_statistical_mechanics()
    
# To directly output the final answer tuple in the required format for submission:
# This reconstruction is to fulfill the specific format requirement.
epsilon_0 = 6.9e-3
T = 4200.0
k_B = 8.617333e-5
N_total_moles = 3.0
num_levels = 5
kbT = k_B * T
n_values = np.arange(1, num_levels + 1)
g_n = 2 * n_values + 1
epsilon_contributions = epsilon_0 * (1 + n_values**n_values)
E_n = np.cumsum(epsilon_contributions)
boltzmann_factors = g_n * np.exp(-E_n / kbT)
Z_1 = np.sum(boltzmann_factors)
moles_in_level = N_total_moles * boltzmann_factors / Z_1
final_answer_tuple = tuple(moles_in_level)
print(f"<<<{final_answer_tuple}>>>")
