import math

def solve_particle_distribution():
    """
    Calculates the partition function and the most probable number of particles 
    in each energy level for a system of distinguishable particles.
    """
    # 1. Define Constants
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T_K = 4200.0  # Kelvin
    k_B_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K
    total_moles = 3.0
    num_levels = 5

    # 2. Calculate Energy Levels (E_n) and Degeneracies (g_n)
    epsilons_eV = []
    energies_eV = []
    degeneracies = []
    
    current_energy_sum_eV = 0.0
    print("--- Calculating Energy Levels and Degeneracies ---")
    for n in range(1, num_levels + 1):
        # Calculate epsilon_n
        epsilon_n = epsilon_0_eV * (1 + n**n)
        epsilons_eV.append(epsilon_n)
        
        # Calculate cumulative energy E_n
        current_energy_sum_eV += epsilon_n
        energies_eV.append(current_energy_sum_eV)
        
        # Calculate degeneracy g_n
        g_n = 2 * n + 1
        degeneracies.append(g_n)
        
        print(f"Level n={n}: E_{n} = {current_energy_sum_eV:.4f} eV, g_{n} = {g_n}")

    # 3. Calculate the Single-Particle Partition Function (z)
    kT_eV = k_B_eV_K * T_K
    
    boltzmann_terms = []
    z = 0.0
    for i in range(num_levels):
        E_n = energies_eV[i]
        g_n = degeneracies[i]
        term = g_n * math.exp(-E_n / kT_eV)
        boltzmann_terms.append(term)
        z += term

    print("\n--- Calculating Partition Function ---")
    print(f"Thermal Energy (kT): {kT_eV:.4f} eV")
    print(f"Single-particle partition function (z): {z:.4f}")
    print("The total partition function Z = z^N, where N = 3 * Avogadro's number.")

    # 4. Calculate the Number of Moles in Each Level
    moles_per_level = []
    for term in boltzmann_terms:
        # Probability P_n = term / z
        # Moles n_n = total_moles * P_n
        moles = total_moles * (term / z)
        moles_per_level.append(moles)
        
    print("\n--- Most Probable Number of Moles per Energy Level ---")
    for i in range(num_levels):
        print(f"Moles in Level E_{i+1}: {moles_per_level[i]:.4f}")

    # 5. Output the final ordered set
    final_answer = tuple(moles_per_level)
    print("\nFinal ordered set of moles (n1, n2, n3, n4, n5):")
    # Using a formatted string to ensure clean tuple output
    formatted_answer = f"({', '.join(f'{x:.4f}' for x in final_answer)})"
    print(formatted_answer)
    
    # The final answer format for the platform
    final_answer_raw = f"<<<{final_answer}>>>"
    # print(final_answer_raw) # This would be for the submission system, not for user display

if __name__ == '__main__':
    solve_particle_distribution()
    # The problem asks for the ordered set (E1, E2, E3, E4, E5) which represents
    # the number of moles.
    # To generate the final output in the required format, we recalculate and format it here.
    epsilon_0_eV = 6.9 * 1e-3
    T_K = 4200.0
    k_B_eV_K = 8.617333262e-5
    total_moles = 3.0
    num_levels = 5
    epsilons_eV = [epsilon_0_eV * (1 + n**n) for n in range(1, num_levels + 1)]
    energies_eV = [sum(epsilons_eV[:i+1]) for i in range(num_levels)]
    degeneracies = [2 * n + 1 for n in range(1, num_levels + 1)]
    kT_eV = k_B_eV_K * T_K
    boltzmann_terms = [g * math.exp(-E / kT_eV) for g, E in zip(degeneracies, energies_eV)]
    z = sum(boltzmann_terms)
    moles_per_level = [total_moles * (term / z) for term in boltzmann_terms]
    
    print(f'<<<{tuple(moles_per_level)}>>>')