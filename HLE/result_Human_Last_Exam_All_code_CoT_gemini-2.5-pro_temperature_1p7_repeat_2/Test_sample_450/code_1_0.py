import numpy as np

def solve_statistical_mechanics():
    """
    Calculates the partition function and the number of moles of particles
    in each of five specified energy levels.
    """
    # --- 1. Define Constants ---
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T_K = 4200.0  # Temperature in Kelvin
    k_B_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K
    total_moles = 3.0 # Total moles of particles (N/N_A)

    # The number of energy levels
    num_levels = 5
    n = np.arange(1, num_levels + 1)

    print("--- Input Parameters ---")
    print(f"Epsilon_0: {epsilon_0_eV:.4e} eV")
    print(f"Temperature T: {T_K} K")
    print(f"Boltzmann constant k_B: {k_B_eV_K:.4e} eV/K")
    k_B_T = k_B_eV_K * T_K
    print(f"Thermal energy k_B*T: {k_B_T:.4f} eV\n")
    
    # --- 2. Calculate Energy Levels and Degeneracies ---
    # Calculate epsilon_n for n=1 to 5
    epsilons = epsilon_0_eV * (1 + n**n)

    # Calculate the energy levels E_n as cumulative sums of epsilon_k
    energies_E_n = np.cumsum(epsilons)
    
    # Calculate the degeneracies g_n = 2n + 1
    degeneracies_g_n = 2 * n + 1
    
    print("--- Calculated Energy Levels and Degeneracies ---")
    for i in range(num_levels):
        print(f"Level n={i+1}: Energy E_{i+1} = {energies_E_n[i]:.4f} eV, Degeneracy g_{i+1} = {degeneracies_g_n[i]}")
    print("")

    # --- 3. Calculate the Single-Particle Partition Function (Z1) ---
    # Calculate Boltzmann factors for each level: g_n * exp(-E_n / k_B*T)
    boltzmann_terms = degeneracies_g_n * np.exp(-energies_E_n / k_B_T)
    
    # The single-particle partition function Z1 is the sum of these terms
    Z1 = np.sum(boltzmann_terms)
    
    print("--- Partition Function ---")
    print("The total partition function Z = (Z1)^N is computationally too large.")
    print(f"The single-particle partition function, Z1 = sum(g_n * exp(-E_n / k_B*T)), is: {Z1:.4f}\n")

    # --- 4. Calculate the Number of Moles in Each Level ---
    # Probability P_n of a particle being in level n
    probabilities_P_n = boltzmann_terms / Z1
    
    # Number of moles n_n in each level
    moles_n_n = total_moles * probabilities_P_n
    
    print("--- Final Result: Moles per Energy Level ---")
    print(f"The most probable number of moles in each energy level (n1, n2, n3, n4, n5) is:")
    
    # Format the final tuple for output
    result_tuple = tuple(moles_n_n)
    print(result_tuple)

    # The final answer in the required format
    print(f"\n<<<({moles_n_n[0]:.4f}, {moles_n_n[1]:.4f}, {moles_n_n[2]:.4f}, {moles_n_n[3]:.4f}, {moles_n_n[4]:.4f})>>>")

if __name__ == '__main__':
    solve_statistical_mechanics()