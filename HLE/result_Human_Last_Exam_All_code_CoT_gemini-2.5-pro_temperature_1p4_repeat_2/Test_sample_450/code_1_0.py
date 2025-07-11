import numpy as np

def solve_particle_distribution():
    """
    Calculates the single-particle partition function and the most probable number of
    particles (in moles) for each of five energy levels in a thermal system.
    """
    # 1. Define Constants
    N_total_moles = 3.0  # Total number of particles in moles (N/N_A)
    epsilon_0 = 6.9e-3      # Base energy in eV
    T = 4200.0              # Temperature in Kelvin
    k_B = 8.617333262e-5    # Boltzmann constant in eV/K

    # Pre-calculate the thermal energy
    kBT = k_B * T

    # 2. Calculate Energy Levels (E_n) and Degeneracies (g_n)
    num_levels = 5
    energies = np.zeros(num_levels)
    degeneracies = np.zeros(num_levels)
    
    current_energy_sum = 0.0
    for n in range(1, num_levels + 1):
        # Calculate the energy contribution for the current level
        epsilon_n = epsilon_0 * (1 + n**n)
        
        # The total energy is a cumulative sum
        current_energy_sum += epsilon_n
        energies[n-1] = current_energy_sum
        
        # Calculate the degeneracy
        degeneracies[n-1] = 2 * n + 1

    # 3. Calculate the Single-Particle Partition Function (q)
    # Z = sum(g_n * exp(-E_n / kBT))
    boltzmann_factors = degeneracies * np.exp(-energies / kBT)
    q = np.sum(boltzmann_factors)

    # 4. Calculate the Probability (P_n) for each Energy Level
    # P_n = (g_n * exp(-E_n / kBT)) / q
    probabilities = boltzmann_factors / q

    # 5. Calculate the Number of Moles in Each Level
    # moles_n = N_total_moles * P_n
    moles_per_level = N_total_moles * probabilities

    # 6. Format the Output
    print(f"The single-particle partition function (q) is: {q}")
    print("The most probable number of particles in each energy level, expressed in moles, is:")
    
    # Create the final ordered set string as requested
    final_result = f"({moles_per_level[0]:.6f}, {moles_per_level[1]:.6f}, {moles_per_level[2]:.6f}, {moles_per_level[3]:.6f}, {moles_per_level[4]:.6f})"
    print(final_result)
    
    # Output the final answer in the specified format
    print(f"<<<{final_result}>>>")

# Execute the function to get the answer
solve_particle_distribution()