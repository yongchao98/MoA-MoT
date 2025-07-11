import math

def solve_statistical_mechanics():
    """
    Solves for the most probable number of particles in moles for each energy level.
    """
    # --- Constants Definition ---
    epsilon0_meV = 6.9  # meV
    epsilon0_eV = epsilon0_meV * 1e-3  # Convert meV to eV
    T = 4200.0  # Temperature in Kelvin
    kB_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K
    N_total_moles = 3.0  # Total number of moles of particles

    # Thermal energy
    kBT = kB_eV_K * T

    # --- Step 1: Calculate Energy Levels (E_n) and Degeneracies (g_n) ---
    num_levels = 5
    energies_E = []
    degeneracies_g = []
    cumulative_energy = 0.0

    print("--- Energy Levels and Degeneracies ---")
    for n in range(1, num_levels + 1):
        # Degeneracy
        g_n = 2 * n + 1
        degeneracies_g.append(g_n)

        # Energy increment
        epsilon_n = epsilon0_eV * (1 + n**n)
        
        # Cumulative Energy E_n
        cumulative_energy += epsilon_n
        energies_E.append(cumulative_energy)

        print(f"Level n={n}: g_{n}={g_n}, E_{n}={cumulative_energy:.4f} eV")

    # --- Step 2: Calculate the Single-Particle Partition Function (Z1) ---
    Z1 = 0.0
    boltzmann_terms = []
    
    print("\n--- Single-Particle Partition Function (Z1) Calculation ---")
    print(f"Z1 = sum(g_n * exp(-E_n / kBT)) where kBT = {kBT:.4f} eV")
    for i in range(num_levels):
        g_n = degeneracies_g[i]
        E_n = energies_E[i]
        term = g_n * math.exp(-E_n / kBT)
        boltzmann_terms.append(term)
        Z1 += term
        print(f"Term for n={i+1}: {g_n} * exp(-{E_n:.4f} / {kBT:.4f}) = {term:.4f}")

    print(f"\nTotal Single-Particle Partition Function Z1 = {Z1:.4f}")

    # --- Step 3: Calculate Moles in Each Level ---
    moles_in_levels = []
    print("\n--- Moles per Energy Level Calculation ---")
    print(f"n_moles_n = (Total Moles) * (g_n * exp(-E_n / kBT)) / Z1")

    for i in range(num_levels):
        n = i + 1
        # Probability P_n = term / Z1
        # Moles = Total Moles * P_n
        moles_n = N_total_moles * (boltzmann_terms[i] / Z1)
        moles_in_levels.append(moles_n)
        
        print(f"Level E_{n}: n_moles = 3.0 * ({boltzmann_terms[i]:.4f} / {Z1:.4f}) = {moles_n:.4f} moles")

    # --- Final Result ---
    final_answer = tuple(moles_in_levels)
    print("\nFinal ordered set of moles in each energy level (E1, E2, E3, E4, E5):")
    print(final_answer)
    print(f"\n<<<{final_answer}>>>")

solve_statistical_mechanics()