import numpy as np

def solve_statistical_mechanics():
    """
    Calculates the partition function and the most probable number of particles
    in each energy level for a given system.
    """
    # 1. Define Constants
    epsilon_0_mev = 6.9
    epsilon_0_ev = epsilon_0_mev * 1e-3  # Convert meV to eV
    T_k = 4200.0  # Temperature in Kelvin
    k_B_ev_k = 8.617333262145e-5  # Boltzmann constant in eV/K
    N_total_moles = 3.0 # Total number of particles is 3 * N_A

    print("--- Constants ---")
    print(f"epsilon_0 = {epsilon_0_ev:.4e} eV")
    print(f"T = {T_k} K")
    print(f"k_B = {k_B_ev_k:.6e} eV/K")
    print(f"Total Moles = {N_total_moles}")
    print("-" * 20)

    # 2. Calculate Energy Levels (E_n) and Degeneracies (g_n)
    n_levels = 5
    epsilons = []
    energies_ev = []
    degeneracies = []
    
    print("--- Energy Levels and Degeneracies ---")
    current_energy = 0.0
    for n in range(1, n_levels + 1):
        # Calculate epsilon_n for the current level
        epsilon_n = epsilon_0_ev * (1 + n**n)
        epsilons.append(epsilon_n)
        
        # Total energy E_n is the cumulative sum of epsilons
        current_energy += epsilon_n
        energies_ev.append(current_energy)
        
        # Calculate degeneracy g_n
        g_n = 2 * n + 1
        degeneracies.append(g_n)
        
        print(f"Level n={n}:")
        print(f"  epsilon_{n} = {epsilon_n:.4e} eV")
        print(f"  Energy E_{n} = {current_energy:.4e} eV")
        print(f"  Degeneracy g_{n} = {g_n}")

    print("-" * 20)

    # 3. Calculate the Single-Particle Partition Function (Z_1)
    beta = 1.0 / (k_B_ev_k * T_k)
    
    print("--- Partition Function Calculation ---")
    print(f"beta = 1 / (k_B * T) = {beta:.4f} eV^-1")

    # Calculate terms for the partition function sum
    z1_terms = [g * np.exp(-beta * E) for g, E in zip(degeneracies, energies_ev)]
    
    # Sum the terms to get the partition function Z_1
    Z_1 = sum(z1_terms)
    
    print(f"The single-particle partition function Z_1 = sum(g_n * exp(-beta * E_n)) is: {Z_1:.4f}")
    print("-" * 20)

    # 4. Calculate the Number of Moles in Each Level
    print("--- Moles per Energy Level ---")
    print("Equation: n_moles = N_total_moles * (g_n * exp(-beta * E_n)) / Z_1")
    
    # The term (g_n * exp(-beta * E_n)) is the same as the z1_terms calculated earlier
    moles_per_level = [(N_total_moles * term) / Z_1 for term in z1_terms]
    
    for i in range(n_levels):
        n = i + 1
        print(f"Moles in level E_{n}: {moles_per_level[i]:.4f}")

    print("-" * 20)
    
    # 5. Output the Final Result
    final_ordered_set = tuple(moles_per_level)
    print("Final ordered set for the number of moles in each energy level (E_1, E_2, E_3, E_4, E_5):")
    print(final_ordered_set)
    
    # The required format for the final answer
    print(f"<<<{final_ordered_set}>>>")

# Execute the function
solve_statistical_mechanics()