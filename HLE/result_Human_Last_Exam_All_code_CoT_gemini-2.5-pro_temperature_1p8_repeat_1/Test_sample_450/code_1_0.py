import math

def solve_statistical_mechanics_problem():
    """
    Solves the statistical mechanics problem as described.
    1. Calculates energy levels and degeneracies.
    2. Computes the single-particle partition function, z.
    3. Determines the most probable number of moles in each energy level.
    """

    # 1. Define Physical Constants
    N_moles = 3.0  # Total number of moles of particles
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T_K = 4200.0  # Temperature in Kelvin
    k_B_eV_K = 8.617333e-5  # Boltzmann constant in eV/K

    # Calculate thermal energy
    k_B_T = k_B_eV_K * T_K

    # 2. Calculate Energy Levels and Degeneracies
    num_levels = 5
    epsilons = []
    energies = []
    degeneracies = []
    
    current_energy_sum = 0.0
    for n in range(1, num_levels + 1):
        # Epsilon_n
        epsilon_n = epsilon_0_eV * (1 + n**n)
        epsilons.append(epsilon_n)
        
        # Cumulative Energy E_n
        current_energy_sum += epsilon_n
        energies.append(current_energy_sum)
        
        # Degeneracy g_n
        g_n = 2 * n + 1
        degeneracies.append(g_n)

    # 3. Calculate the Single-Particle Partition Function (z)
    z = 0.0
    boltzmann_terms = []
    for i in range(num_levels):
        g_i = degeneracies[i]
        E_i = energies[i]
        term = g_i * math.exp(-E_i / k_B_T)
        boltzmann_terms.append(term)
        z += term
        
    print(f"The single-particle partition function z is: {z:.4f}")
    # Note: The total partition function Z = z^N is computationally infeasible to calculate.
    
    # 4. Calculate the Number of Moles in Each Energy Level (n_i)
    moles_per_level = []
    for i in range(num_levels):
        n_i = N_moles * boltzmann_terms[i] / z
        moles_per_level.append(n_i)
    
    # 5. Output the Results
    print("\nThe most probable number of moles in each energy level is:")
    print(f"(n1, n2, n3, n4, n5) = ({moles_per_level[0]:.4f}, {moles_per_level[1]:.4f}, {moles_per_level[2]:.4f}, {moles_per_level[3]:.4f}, {moles_per_level[4]:.4f})")
    
    # Format final answer for direct extraction
    final_answer_tuple = tuple(moles_per_level)
    print(f"\n<<<({final_answer_tuple[0]:.4f}, {final_answer_tuple[1]:.4f}, {final_answer_tuple[2]:.4f}, {final_answer_tuple[3]:.4f}, {final_answer_tuple[4]:.4f})>>>")

# Execute the function
solve_statistical_mechanics_problem()