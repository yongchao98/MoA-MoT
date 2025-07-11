import numpy as np

def solve_statistical_mechanics():
    """
    Calculates the partition function and the number of moles in each energy level
    for a system of distinguishable particles.
    """
    # 1. Define Constants and Parameters
    N_total_moles = 3.0
    epsilon_0_meV = 6.9
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T_K = 4200.0  # Temperature in Kelvin
    k_B_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K

    # Calculate thermal energy
    kBT_eV = k_B_eV_K * T_K

    print("--- System Parameters ---")
    print(f"Total number of moles: {N_total_moles}")
    print(f"Epsilon_0: {epsilon_0_eV:.5f} eV")
    print(f"Temperature T: {T_K} K")
    print(f"Thermal Energy (k_B * T): {kBT_eV:.5f} eV\n")

    # 2. Calculate Energy Levels and Degeneracies
    n_values = np.arange(1, 6)
    
    # Degeneracies: g_n = 2n + 1
    g_n = 2 * n_values + 1
    
    # Energy increments: epsilon_n = epsilon_0 * (1 + n^n)
    delta_epsilon_n = epsilon_0_eV * (1 + n_values**n_values)
    
    # Cumulative energy levels E_n
    E_n = np.cumsum(delta_epsilon_n)

    print("--- Energy Levels and Degeneracies ---")
    print("Level (n) | Degeneracy (g_n) | Energy (E_n, eV)")
    print("---------------------------------------------")
    for i in range(5):
        print(f"    {n_values[i]:<5} | {g_n[i]:<16} | {E_n[i]:<18.5f}")
    print("\n")

    # 3. Calculate Partition Function Z
    # Z = sum(g_n * exp(-E_n / kBT))
    boltzmann_terms = g_n * np.exp(-E_n / kBT_eV)
    Z = np.sum(boltzmann_terms)

    print("--- Partition Function (Z) Calculation ---")
    print("Z = g_1*exp(-E_1/kBT) + g_2*exp(-E_2/kBT) + ... + g_5*exp(-E_5/kBT)")
    equation_str = "Z = " + " + ".join([f"{term:.4f}" for term in boltzmann_terms])
    print(equation_str)
    print(f"Z = {Z:.4f}\n")

    # 4. Calculate Number of Moles in Each Level
    # n_i = N_total_moles * (g_i * exp(-E_i / kBT)) / Z
    moles_in_level = N_total_moles * boltzmann_terms / Z

    print("--- Moles per Energy Level (n_i) Calculation ---")
    print("Formula: n_i = (Total Moles) * (g_i * exp(-E_i / kBT)) / Z")
    for i in range(5):
        print(f"n_{i+1} = {N_total_moles:.1f} * ({boltzmann_terms[i]:.4f}) / {Z:.4f} = {moles_in_level[i]:.4f} moles")

    # Final ordered set
    result_tuple = tuple(moles_in_level)
    print("\nFinal ordered set for the number of moles in each energy level (E_1, E_2, E_3, E_4, E_5):")
    print(f"({result_tuple[0]:.4f}, {result_tuple[1]:.4f}, {result_tuple[2]:.4f}, {result_tuple[3]:.4f}, {result_tuple[4]:.4f})")
    
    # Format the final answer string as requested
    final_answer_str = f"({result_tuple[0]:.4f}, {result_tuple[1]:.4f}, {result_tuple[2]:.4f}, {result_tuple[3]:.4f}, {result_tuple[4]:.4f})"
    return final_answer_str

# Execute the function and print the final answer in the required format
final_answer = solve_statistical_mechanics()
print(f"\n<<<{final_answer}>>>")