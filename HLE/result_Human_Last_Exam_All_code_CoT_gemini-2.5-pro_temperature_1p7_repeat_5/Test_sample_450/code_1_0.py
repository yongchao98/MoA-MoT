import math

def solve_statistical_mechanics():
    """
    Calculates the single-particle partition function and the most probable number of moles
    in each energy level for a system of distinguishable particles.
    """
    # 1. Define Constants and Parameters
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T = 4200.0  # Temperature in Kelvin
    kB_eV_K = 8.617333e-5  # Boltzmann constant in eV/K
    N_moles_total = 3.0  # Total number of moles

    # Calculate thermal energy for the exponent
    kT_eV = kB_eV_K * T

    # 2. Calculate Energy Levels and Degeneracies
    n_values = range(1, 6)
    
    # Degeneracies g_n = 2n + 1
    g_n = [2 * n + 1 for n in n_values]

    # Epsilon components epsilon_n = epsilon_0 * (1 + n^n)
    epsilon_n_components = [epsilon_0_eV * (1 + n**n) for n in n_values]

    # Cumulative energy levels E_n = sum of epsilon components
    E_n = []
    current_energy = 0
    for eps in epsilon_n_components:
        current_energy += eps
        E_n.append(current_energy)

    # 3. Calculate the Single-Particle Partition Function (Z1)
    boltzmann_terms = []
    for i in range(len(n_values)):
        exponent = -E_n[i] / kT_eV
        term = g_n[i] * math.exp(exponent)
        boltzmann_terms.append(term)
    
    Z1 = sum(boltzmann_terms)

    # 4. Calculate the Number of Moles per Level
    probabilities = [term / Z1 for term in boltzmann_terms]
    moles_per_level = [N_moles_total * p for p in probabilities]

    # 5. Output the results
    print("--- System Parameters and Constants ---")
    print(f"Epsilon_0: {epsilon_0_eV:.4f} eV")
    print(f"Temperature (T): {T} K")
    print(f"Thermal Energy (k_B * T): {kT_eV:.4f} eV")
    
    print("\n--- Energy Levels and Partition Function ---")
    print("Level (n) | Degeneracy (g_n) | Energy (E_n) [eV]   | Partition Function Term")
    print("----------|------------------|--------------------|-------------------------")
    for i, n in enumerate(n_values):
        print(f"    {n}     | {g_n[i]:<16} | {E_n[i]:<18.4f} | {boltzmann_terms[i]:.4f}")
    
    print(f"\nTotal Single-Particle Partition Function (Z1): {Z1:.4f}")

    print("\n--- Final Result: Moles per Energy Level ---")
    print("The final equation for moles in level n is: n_n = N_total_moles * (g_n * exp(-E_n / (k_B * T))) / Z1")
    print("\nExample for Level 1:")
    print(f"n_1 = {N_moles_total:.1f} * ({g_n[0]} * exp(-{E_n[0]:.4f} / {kT_eV:.4f})) / {Z1:.4f} = {moles_per_level[0]:.4f} moles")
    
    print("\nOrdered set of moles for levels (E1, E2, E3, E4, E5):")
    final_tuple = tuple(moles_per_level)
    print(final_tuple)
    
    # Required final answer format
    print(f"\n<<<{final_tuple}>>>")

solve_statistical_mechanics()