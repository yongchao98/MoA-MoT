import math

def solve_statistical_mechanics_problem():
    """
    Solves the statistical mechanics problem to find the number of moles in each energy level.
    """
    # Step 1: Define Constants
    epsilon_0_meV = 6.9
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T_K = 4200.0  # Temperature in Kelvin
    k_B_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K
    N_total_moles = 3.0 # Total number of moles of particles
    num_levels = 5

    # Thermal energy in eV
    kBT_eV = k_B_eV_K * T_K

    # Step 2: Calculate Energy Levels and Degeneracies
    epsilons = []
    energies_cumulative = []
    degeneracies = []
    current_energy_sum = 0.0

    print("--- Calculating Energy Levels and Degeneracies ---")
    print(f"{'Level n':<10}{'epsilon_n (eV)':<20}{'E_n (eV)':<20}{'g_n':<10}")
    print("-" * 65)

    for n in range(1, num_levels + 1):
        # Epsilon for this level
        epsilon_n = epsilon_0_eV * (1 + n**n)
        
        # Cumulative energy E_n
        current_energy_sum += epsilon_n
        
        # Degeneracy g_n
        g_n = 2 * n + 1
        
        energies_cumulative.append(current_energy_sum)
        degeneracies.append(g_n)
        
        print(f"{n:<10}{epsilon_n:<20.6e}{current_energy_sum:<20.6e}{g_n:<10}")
    print("\n")

    # Step 3: Calculate Single-Particle Partition Function Z1
    boltzmann_terms = []
    Z1 = 0.0

    for i in range(num_levels):
        E_n = energies_cumulative[i]
        g_n = degeneracies[i]
        exponent = -E_n / kBT_eV
        term = g_n * math.exp(exponent)
        boltzmann_terms.append(term)
        Z1 += term

    print("--- Calculating Single-Particle Partition Function (Z1) ---")
    print(f"Using Z1 = sum(g_n * exp(-E_n / kBT)), where kBT = {kBT_eV:.6f} eV")
    print("-" * 65)
    for i in range(num_levels):
        print(f"Term for n={i+1}: {degeneracies[i]} * exp(-{energies_cumulative[i]:.6e} / {kBT_eV:.6f}) = {boltzmann_terms[i]:.6f}")
    print("-" * 65)
    print(f"Total Single-Particle Partition Function Z1 = {Z1:.6f}\n")


    # Step 4: Calculate the Number of Moles in Each Level
    moles_per_level = []
    print("--- Calculating Moles per Energy Level ---")
    print(f"Using n_moles = N_total_moles * (Term / Z1)")
    print(f"Total Moles = {N_total_moles}")
    print("-" * 65)
    
    for i in range(num_levels):
        n_moles = N_total_moles * boltzmann_terms[i] / Z1
        moles_per_level.append(n_moles)
        # Outputting the numbers in the final equation as requested
        print(f"Moles in E_{i+1} = {N_total_moles:.1f} * ({boltzmann_terms[i]:.6f} / {Z1:.6f}) = {n_moles:.6f}")
    
    print("-" * 65)
    print("\n")

    # Step 5: Provide the Final Answer
    # The problem requests the ordered set (E1, E2, E3, E4, E5) representing the number of moles.
    final_answer_tuple = tuple(moles_per_level)
    print("Final Answer: The ordered set of the number of moles in each energy level (E1 to E5) is:")
    print(final_answer_tuple)
    print(f"\n<<<{final_answer_tuple}>>>")

if __name__ == '__main__':
    solve_statistical_mechanics_problem()