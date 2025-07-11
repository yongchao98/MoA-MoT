import numpy as np

def solve_statistical_mechanics():
    """
    Calculates the partition function and the most probable number of moles
    in each energy level for a given physical system.
    """
    # 1. Define Constants
    epsilon_0 = 6.9e-3  # eV
    T = 4200.0          # K
    k_B = 8.617333e-5   # eV/K
    N_moles = 3.0       # Total number of moles
    n_levels = 5

    # 2. Calculate Energy Levels and Degeneracies
    energies = []
    degeneracies = []
    cumulative_epsilon = 0.0
    for n in range(1, n_levels + 1):
        # Epsilon_n calculation
        epsilon_n = epsilon_0 * (1 + n**n)
        # E_n calculation (cumulative sum)
        cumulative_epsilon += epsilon_n
        energies.append(cumulative_epsilon)
        # g_n calculation
        g_n = 2 * n + 1
        degeneracies.append(g_n)

    # 3. Calculate the Partition Function (Z_1)
    kBT = k_B * T
    Z_1 = 0.0
    boltzmann_terms = []
    for i in range(n_levels):
        exponent = -energies[i] / kBT
        term = degeneracies[i] * np.exp(exponent)
        boltzmann_terms.append(term)
        Z_1 += term

    print("The single-particle partition function Z is calculated as:")
    z_string_values = "Z = "
    for i in range(n_levels):
        z_string_values += f"{degeneracies[i]} * exp(-{energies[i]:.4f} / {kBT:.4f})"
        if i < n_levels - 1:
            z_string_values += " + "
    print(z_string_values)
    
    z_sum_string = "Z = "
    for i, term in enumerate(boltzmann_terms):
        z_sum_string += f"{term:.4f}"
        if i < len(boltzmann_terms) - 1:
            z_sum_string += " + "
    print(z_sum_string)
    print(f"Z = {Z_1:.4f}\n")
    

    # 4. Calculate Occupation Numbers in Moles
    occupation_moles = []
    print("The most probable number of moles in each level n_i is calculated as:")
    for i in range(n_levels):
        n_i_moles = N_moles * boltzmann_terms[i] / Z_1
        occupation_moles.append(n_i_moles)
        print(f"n_{i+1} = {N_moles} * ({boltzmann_terms[i]:.4f}) / {Z_1:.4f} = {n_i_moles:.4f} moles")

    # 5. Output Final Answer
    print("\nThe ordered set representing the number of moles in each energy level (E_1, E_2, E_3, E_4, E_5) is:")
    final_tuple_str = ", ".join([f"{x:.4f}" for x in occupation_moles])
    print(f"({final_tuple_str})")
    
    # Final answer in the required format
    final_answer = tuple(occupation_moles)
    print(f"\n<<<{final_answer}>>>")

solve_statistical_mechanics()