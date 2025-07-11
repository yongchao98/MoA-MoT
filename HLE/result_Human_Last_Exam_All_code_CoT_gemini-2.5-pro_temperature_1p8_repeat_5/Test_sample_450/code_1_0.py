import math

def solve_statistical_mechanics():
    """
    Calculates the single-particle partition function and the molar distribution
    of particles across five energy levels.
    """
    # 1. Define Constants
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert to eV
    T_K = 4200  # Kelvin
    k_B_eV_K = 8.617333262145e-5  # Boltzmann constant in eV/K
    N_total_moles = 3.0  # Total number of moles

    # Calculate thermal energy
    kBT = k_B_eV_K * T_K

    # 2. & 3. Calculate Energy Levels (E_n) and Degeneracies (g_n)
    energies = []
    degeneracies = []
    cumulative_energy_eV = 0.0

    for n in range(1, 6):
        # Calculate specific energy epsilon_n
        eps_n = epsilon_0_eV * (1 + n**n)
        
        # Calculate cumulative energy E_n
        cumulative_energy_eV += eps_n
        energies.append(cumulative_energy_eV)
        
        # Calculate degeneracy g_n
        g_n = 2 * n + 1
        degeneracies.append(g_n)

    # 4. Calculate Single-Particle Partition Function (Z_1)
    Z1 = 0.0
    boltzmann_terms = []
    for i in range(5):
        E_n = energies[i]
        g_n = degeneracies[i]
        term = g_n * math.exp(-E_n / kBT)
        boltzmann_terms.append(term)
        Z1 += term
        
    print(f"The single-particle partition function Z_1 is: {Z1}")

    # 5. Calculate the most probable number of moles in each energy level
    moles_distribution = []
    for term in boltzmann_terms:
        # Probability P_n = term / Z1
        # Moles_n = N_total_moles * P_n
        moles_n = N_total_moles * (term / Z1)
        moles_distribution.append(moles_n)
        
    # 6. Output the final result
    print("The most probable number of particles in each energy level, in moles, is:")
    print(f"(E1, E2, E3, E4, E5) = {tuple(moles_distribution)}")

solve_statistical_mechanics()
<<< (0.795604118744093, 1.2054238555819777, 0.9897282860829895, 0.00924160352226871, 2.136068671158586e-26) >>>