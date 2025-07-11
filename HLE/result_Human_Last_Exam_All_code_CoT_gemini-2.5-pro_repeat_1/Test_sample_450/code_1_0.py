import math

def solve_statistical_mechanics():
    """
    Calculates the partition function and the most probable number of moles
    of particles in each of five energy levels.
    """
    # 1. Define Constants
    epsilon_0_meV = 6.9  # meV
    T = 4200.0  # K
    total_moles = 3.0

    # Physical constants
    eV_to_J = 1.60218e-19  # J/eV
    k_B_J_K = 1.380649e-23  # J/K
    
    # Convert epsilon_0 to Joules
    epsilon_0 = epsilon_0_meV * 1e-3 * eV_to_J

    # Thermal energy
    kT = k_B_J_K * T

    # 2. Calculate Energy Levels and Degeneracies
    epsilons = []
    energies = []
    degeneracies = []
    current_energy = 0.0

    for n in range(1, 6):
        # Epsilon_n for the current level
        epsilon_n = epsilon_0 * (1 + n**n)
        epsilons.append(epsilon_n)
        
        # Cumulative energy E_n
        current_energy += epsilon_n
        energies.append(current_energy)
        
        # Degeneracy g_n
        g_n = 2 * n + 1
        degeneracies.append(g_n)

    # 3. Calculate the Single-Particle Partition Function (Z1)
    Z1 = 0.0
    boltzmann_terms = []
    for i in range(5):
        g_n = degeneracies[i]
        E_n = energies[i]
        term = g_n * math.exp(-E_n / kT)
        boltzmann_terms.append(term)
        Z1 += term
        
    print(f"The single-particle partition function Z_1 is: {Z1:.4f}")
    print("-" * 20)

    # 4. Calculate the Number of Moles per Level
    moles_per_level = []
    print("Calculation of the number of moles in each energy level:")
    for i in range(5):
        n_level = i + 1
        g_n = degeneracies[i]
        E_n = energies[i]
        
        # Probability of a particle being in level n
        probability = boltzmann_terms[i] / Z1
        
        # Most probable number of moles in level n
        moles_n = total_moles * probability
        moles_per_level.append(moles_n)
        
        # Output the equation with numbers
        print(f"\nMoles in level E_{n_level}:")
        print(f"n_{n_level} = (Total Moles) * g_{n_level} * exp(-E_{n_level} / kT) / Z_1")
        print(f"n_{n_level} = {total_moles:.2f} * {g_n} * exp(-{E_n:.4e} / {kT:.4e}) / {Z1:.4f}")
        print(f"n_{n_level} = {moles_n:.6f} moles")

    # 5. Output the final ordered set
    final_answer = tuple(moles_per_level)
    print("\n" + "-" * 20)
    print("The ordered set of the number of moles in each energy level (E1, E2, E3, E4, E5) is:")
    print(final_answer)
    
    # Final answer in the required format
    print(f"\n<<<{final_answer}>>>")


solve_statistical_mechanics()