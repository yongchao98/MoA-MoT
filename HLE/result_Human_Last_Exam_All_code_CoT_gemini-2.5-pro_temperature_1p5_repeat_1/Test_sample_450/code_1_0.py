import numpy as np

def solve_statistical_mechanics():
    """
    Calculates the single-particle partition function and the most probable number of
    moles in each of the five energy levels for a system of distinguishable particles.
    """
    # 1. Define Constants
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert to eV
    T_K = 4200  # Kelvin
    k_B_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K
    total_moles = 3.0

    # Calculate kT product
    kT_eV = k_B_eV_K * T_K

    # 2. Calculate Energy Levels and Degeneracies
    n_values = np.arange(1, 6)
    degeneracies_g = 2 * n_values + 1
    
    epsilons_n = epsilon_0_eV * (1 + n_values**n_values)
    energies_E = np.cumsum(epsilons_n)

    print("--- System Parameters ---")
    print(f"epsilon_0 = {epsilon_0_eV:.5f} eV")
    print(f"Temperature T = {T_K} K")
    print(f"k_B * T = {kT_eV:.5f} eV")
    print("\n--- Energy Levels and Degeneracies ---")
    for i, n in enumerate(n_values):
        print(f"Level n={n}: E_{n} = {energies_E[i]:.5f} eV, g_{n} = {degeneracies_g[i]}")

    # 3. Calculate the Single-Particle Partition Function (Z_1)
    boltzmann_factors = np.exp(-energies_E / kT_eV)
    terms_z1 = degeneracies_g * boltzmann_factors
    z1 = np.sum(terms_z1)

    print("\n--- Calculation of Single-Particle Partition Function Z_1 ---")
    print("Z_1 = sum(g_n * exp(-E_n / (k_B * T)))")
    for i, n in enumerate(n_values):
        print(f"Term for n={n}: {degeneracies_g[i]} * exp(-{energies_E[i]:.5f} / {kT_eV:.5f}) = {terms_z1[i]:.5f}")
    print(f"\nSingle-Particle Partition Function Z_1 = {z1:.5f}")

    # 4. Calculate the Number of Moles in Each Energy Level
    probabilities = terms_z1 / z1
    moles_in_levels = total_moles * probabilities

    print("\n--- Calculation of Moles in Each Energy Level ---")
    print("moles_n = Total Moles * (g_n * exp(-E_n / (k_B * T))) / Z_1")
    for i, n in enumerate(n_values):
        print(f"Moles in E_{n}: {total_moles:.1f} * ({degeneracies_g[i]} * exp(-{energies_E[i]:.5f} / {kT_eV:.5f})) / {z1:.5f} = {moles_in_levels[i]:.5f} moles")
    
    print(f"\nCheck: Total moles = {np.sum(moles_in_levels):.5f}")

    # Final answer formatting
    final_answer = tuple(moles_in_levels)
    print("\nFinal ordered set of moles in each energy level (E1, E2, E3, E4, E5):")
    print(f"<<<{final_answer}>>>")


solve_statistical_mechanics()