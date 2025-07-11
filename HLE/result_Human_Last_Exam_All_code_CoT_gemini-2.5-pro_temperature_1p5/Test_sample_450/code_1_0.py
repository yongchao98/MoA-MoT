import math

def solve_statistical_mechanics_problem():
    """
    Solves the statistical mechanics problem as described.
    Calculates the single-particle partition function and the number of moles
    of particles in each of the five energy levels.
    """
    # Step 1: Define Constants
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV * 1e-3  # convert to eV
    T_K = 4200.0  # Kelvin
    k_B_eV_K = 8.617333e-5  # Boltzmann constant in eV/K
    total_moles = 3.0

    print("--- Constants and Parameters ---")
    print(f"epsilon_0 = {epsilon_0_eV:.6f} eV")
    print(f"Temperature T = {T_K} K")
    print(f"Boltzmann constant k_B = {k_B_eV_K} eV/K")
    print(f"Total number of moles = {total_moles}")
    print("-" * 30)

    # Calculate thermal energy kT
    kT = k_B_eV_K * T_K
    print(f"Thermal energy kT = {kT:.6f} eV\n")

    # Step 2: Calculate Energy Levels and Degeneracies
    energies_E = []
    degeneracies_g = []
    epsilons = []
    
    current_E = 0.0
    for n in range(1, 6):
        # Calculate degeneracy
        g_n = 2 * n + 1
        degeneracies_g.append(g_n)

        # Calculate epsilon_n
        epsilon_n = epsilon_0_eV * (1 + n**n)
        epsilons.append(epsilon_n)
        
        # Calculate cumulative energy E_n
        current_E += epsilon_n
        energies_E.append(current_E)

    # Step 3: Compute the Single-Particle Partition Function (Z_1)
    Z_single = 0.0
    boltzmann_terms = []
    
    print("--- Calculation for Partition Function Z_1 ---")
    print("Z_1 = sum(g_n * exp(-E_n / kT)) for n=1 to 5")
    print("-" * 45)
    print(f"{'Level (n)':<10} {'g_n':<5} {'E_n (eV)':<12} {'g_n * exp(-E_n/kT)':<20}")
    print("-" * 45)

    for i in range(5):
        n = i + 1
        E_n = energies_E[i]
        g_n = degeneracies_g[i]
        
        # Calculate Boltzmann factor term: g_n * exp(-E_n / kT)
        term = g_n * math.exp(-E_n / kT)
        boltzmann_terms.append(term)
        Z_single += term
        
        print(f"{n:<10} {g_n:<5} {E_n:<12.6f} {term:<20.6f}")

    print("-" * 45)
    print(f"The single-particle partition function Z_1 = {Z_single:.6f}\n")

    # Step 4: Calculate Occupation Numbers in Moles
    moles_per_level = []
    
    print("--- Calculation for Number of Moles per Level ---")
    print("n_moles = Total Moles * (g_n * exp(-E_n / kT)) / Z_1")
    print("-" * 45)
    print(f"{'Level (n)':<10} {'n_moles':<15}")
    print("-" * 45)

    for i in range(5):
        n = i + 1
        # n_moles = total_moles * P_n = total_moles * (term / Z_single)
        n_moles = total_moles * boltzmann_terms[i] / Z_single
        moles_per_level.append(n_moles)
        print(f"{n:<10} {n_moles:<15.6f}")

    print("-" * 45)

    # Step 5: Format the Output
    final_moles_tuple = tuple(moles_per_level)
    print("\nFinal Result:")
    print("The most probable number of particles in each energy level, in moles, is:")
    print(f"(moles_E1, moles_E2, moles_E3, moles_E4, moles_E5) = {final_moles_tuple}")
    
    # Final answer in the required format
    print(f"\n<<<{final_moles_tuple}>>>")

if __name__ == '__main__':
    solve_statistical_mechanics_problem()