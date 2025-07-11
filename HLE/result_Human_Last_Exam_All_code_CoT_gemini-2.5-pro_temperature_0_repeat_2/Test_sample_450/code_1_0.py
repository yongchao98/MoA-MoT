import math

def solve_statistical_mechanics():
    """
    Calculates the single-particle partition function and the number of moles
    of particles in each of five energy levels for a given system.
    """
    # Step 1: Define constants
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T_K = 4200.0  # Kelvin
    k_B_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K
    N_total_moles = 3.0 # Total number of particles is 3 * N_A

    # Thermal energy
    k_B_T = k_B_eV_K * T_K

    # Step 2: Calculate energy levels and degeneracies
    energies = []
    degeneracies = []
    current_energy_E = 0.0
    
    print("Calculating Energy Levels and Degeneracies:")
    for n in range(1, 6):
        epsilon_n = epsilon_0_eV * (1 + n**n)
        current_energy_E += epsilon_n
        g_n = 2 * n + 1
        
        energies.append(current_energy_E)
        degeneracies.append(g_n)
        
        print(f"Level n={n}: E_{n} = {current_energy_E:.4f} eV, g_{n} = {g_n}")

    # Step 3: Calculate the single-particle partition function (Z_1)
    Z1 = 0.0
    partition_terms = []
    for i in range(5):
        E_n = energies[i]
        g_n = degeneracies[i]
        term = g_n * math.exp(-E_n / k_B_T)
        partition_terms.append(term)
        Z1 += term

    print(f"\nThe single-particle partition function Z_1 = {Z1:.4f}")

    # Step 4: Calculate the most probable number of moles in each level
    moles_per_level = []
    for term in partition_terms:
        moles_n = N_total_moles * term / Z1
        moles_per_level.append(moles_n)

    print("\nMost probable number of particles in each energy level (in moles):")
    for i in range(5):
        print(f"Level E_{i+1}: {moles_per_level[i]:.4f} moles")
        
    # Step 5: Format the final answer
    final_answer = tuple(moles_per_level)
    print(f"\nFinal ordered set of moles (E1, E2, E3, E4, E5):")
    print(final_answer)
    
    # The final answer in the required format
    # Using f-string formatting to control precision in the final output string
    answer_str = f"({', '.join(f'{x:.4f}' for x in final_answer)})"
    print(f"\n<<<{answer_str}>>>")

solve_statistical_mechanics()