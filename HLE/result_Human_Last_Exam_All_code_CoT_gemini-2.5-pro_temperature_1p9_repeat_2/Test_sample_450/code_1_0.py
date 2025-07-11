import math

def solve_statistical_mechanics_problem():
    """
    Solves for the partition function and the most probable number of particles
    in each energy level for a given statistical mechanics system.
    """

    # 1. Define Constants and Parameters
    epsilon_0_meV = 6.9
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T_K = 4200.0  # Temperature in Kelvin
    k_B_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K
    N_total_moles = 3.0

    print("--- Constants and Parameters ---")
    print(f"epsilon_0 = {epsilon_0_eV:.3e} eV")
    print(f"Temperature T = {T_K} K")
    print(f"Boltzmann constant k_B = {k_B_eV_K:.6e} eV/K")
    print(f"Total number of moles = {N_total_moles}")
    print("-" * 30)

    # Calculate thermal energy kT for convenience
    kT = k_B_eV_K * T_K
    print(f"Thermal energy kT = {kT:.4f} eV\n")

    # 2. Calculate Energy Levels and Degeneracies
    energy_levels = []
    degeneracies = []
    cumulative_energy_term = 0.0

    print("--- Energy Levels (En) and Degeneracies (gn) ---")
    for n in range(1, 6):
        # Calculate the energy increment for level n
        epsilon_n = epsilon_0_eV * (1 + n**n)
        
        # The energy levels E_n are cumulative sums
        cumulative_energy_term += epsilon_n
        energy_levels.append(cumulative_energy_term)
        
        # Calculate degeneracy for level n
        g_n = 2 * n + 1
        degeneracies.append(g_n)

        print(f"For n={n}: E_{n} = {cumulative_energy_term:.4e} eV, g_{n} = {g_n}")
    print("-" * 30)
    
    # 3. Calculate the Single-Particle Partition Function (z)
    z = 0.0
    z_terms = []
    print("--- Single-Particle Partition Function (z) ---")
    print("z = sum[ gn * exp(-En / kT) ]")
    
    equation_str = "z = "
    for i in range(len(energy_levels)):
        g = degeneracies[i]
        E = energy_levels[i]
        
        boltzmann_factor = math.exp(-E / kT)
        term = g * boltzmann_factor
        z_terms.append(term)
        z += term
        
        equation_str += f"{term:.4f}"
        if i < len(energy_levels) - 1:
            equation_str += " + "
    
    print(f"The calculation is:\n{equation_str}")
    print(f"The single-particle partition function z = {z:.4f}")
    print("The total partition function is Z = z^N for N distinguishable particles.")
    print("-" * 30)
    
    # 4. Calculate the Number of Moles in Each Level
    moles_in_levels = []
    print("--- Most Probable Number of Moles per Level ---")
    print("moles_n = N_total_moles * (gn * exp(-En / kT)) / z")

    for i in range(len(z_terms)):
        moles_n = N_total_moles * z_terms[i] / z
        moles_in_levels.append(moles_n)
        print(f"moles_{i+1} = {N_total_moles} * ({z_terms[i]:.4f} / {z:.4f}) = {moles_n:.4f} moles")
    print("-" * 30)

    # 5. Format and Output Final Answer
    final_answer_tuple_str = f"({moles_in_levels[0]:.4f}, {moles_in_levels[1]:.4f}, {moles_in_levels[2]:.4f}, {moles_in_levels[3]:.4f}, {moles_in_levels[4]:.4f})"
    print(f"The final ordered set of moles in each energy level (E1, E2, E3, E4, E5) is:")
    print(final_answer_tuple_str)

    # Final answer in the required format
    final_answer_str = f"({', '.join(f'{m:.4f}' for m in moles_in_levels)})"
    print(f"\n<<<{final_answer_str}>>>")

solve_statistical_mechanics_problem()