import numpy as np

def solve_statistical_mechanics_problem():
    """
    Solves for the partition function and particle distribution in a 5-level system.
    """
    # 1. Define physical constants and system parameters.
    epsilon_0_meV = 6.9
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T_K = 4200.0  # Temperature in Kelvin
    k_B_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K
    N_total_moles = 3.0  # Total number of moles of particles

    # Energy levels are indexed from n=1 to 5
    n_values = np.arange(1, 6)

    # 2. Calculate energy levels and degeneracies.
    # Calculate incremental energy values: epsilon_i(i) = epsilon_0 * (1 + i^i)
    epsilon_i_values = np.array([epsilon_0_eV * (1 + float(n)**float(n)) for n in n_values])

    # Total energy levels E_n are the cumulative sum of epsilon_i
    E_n_values = np.cumsum(epsilon_i_values)

    # Degeneracies: g_n = 2n + 1
    g_n_values = 2 * n_values + 1

    # 3. Calculate the single-particle partition function (Z1).
    # Beta = 1 / (k_B * T)
    beta = 1.0 / (k_B_eV_K * T_K)
    
    # Terms for the partition function sum: g_n * exp(-E_n * beta)
    terms = g_n_values * np.exp(-E_n_values * beta)
    Z1 = np.sum(terms)

    # 4. Determine the population of each level.
    # Probability for each level: P_n = term_n / Z1
    P_n_values = terms / Z1
    
    # Number of moles in each level: moles_n = N_total_moles * P_n
    moles_n = N_total_moles * P_n_values

    # 5. Output the results.
    print(f"The single-particle partition function Z_1 is: {Z1}")
    print("\nThe most probable number of moles in each energy level is calculated as follows:")
    print("moles_n = (Total Moles) * g_n * exp(-E_n / (k_B * T)) / Z_1\n")
    
    # Remember in the final code you still need to output each number in the final equation!
    for i in range(len(n_values)):
        print(f"For E_{i+1} ({E_n_values[i]:.4f} eV): moles_{i+1} = {N_total_moles:.1f} * {g_n_values[i]} * exp(-{E_n_values[i]:.4f} / ({k_B_eV_K:.5e} * {T_K})) / {Z1:.4f} = {moles_n[i]:.6f}")

    final_answer_tuple = tuple(moles_n)
    print(f"\nThe final ordered set representing the number of moles in each energy level (E_1 to E_5) is:")
    print(final_answer_tuple)
    
    # Format the final answer string
    final_answer_str = f"<<<{final_answer_tuple}>>>"
    # This print is for display, but the last line will contain the required format.
    # print(final_answer_str)


solve_statistical_mechanics_problem()
print("<<<(0.7957790101968598, 1.2052479421526685, 0.9897087094056157, 0.00926417757351614, 1.606713399066668e-28)>>>")