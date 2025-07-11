import numpy as np

def solve_particle_distribution():
    """
    Calculates the partition function and the most probable number of particles 
    in each energy level for a given system.
    """
    # 1. Define Physical Constants
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert to eV
    T = 4200.0  # Temperature in Kelvin
    k_B_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K
    n_total_moles = 3.0 # Total number of moles

    # Product of k_B and T, the thermal energy
    kT = k_B_eV_K * T

    # Create an array for the level indices n = 1, 2, 3, 4, 5
    n_levels = np.arange(1, 6)

    # 2. Calculate Energy Levels (E_n)
    # Calculate individual energy contributions epsilon_n
    epsilon_n = epsilon_0_eV * (1 + n_levels**n_levels)
    # Calculate the cumulative energy levels E_n
    E_n = np.cumsum(epsilon_n)

    # 3. Calculate Degeneracies (g_n)
    g_n = 2 * n_levels + 1

    # 4. Calculate the Single-Particle Partition Function (Z_1)
    # Calculate the terms in the partition function sum
    boltzmann_terms = g_n * np.exp(-E_n / kT)
    # Sum the terms to get the partition function Z_1
    Z1 = np.sum(boltzmann_terms)

    # 5. Calculate the Number of Moles in Each Level (n_n)
    # The probability of a particle being in level n is P_n = boltzmann_terms[n-1] / Z1
    # The number of moles in level n is n_n = n_total_moles * P_n
    moles_per_level = n_total_moles * boltzmann_terms / Z1

    # 6. Format and Print the Output
    print(f"The thermal energy k_B * T = {kT:.5f} eV")
    print(f"The single-particle partition function Z_1 = {Z1:.5f}\n")
    print("The total number of moles is 3.0.")
    print("The number of moles in each energy level is calculated as:")
    print("n_i = n_total * (g_i * exp(-E_i / (k_B * T))) / Z_1\n")
    
    for i in range(len(n_levels)):
        n = n_levels[i]
        print(f"--- Level n={n} ---")
        print(f"Energy E_{n} = {E_n[i]:.5f} eV")
        print(f"Degeneracy g_{n} = {g_n[i]}")
        # Final equation components
        print(f"Calculation: n_{n} = {n_total_moles} * ({g_n[i]} * exp(-{E_n[i]:.5f} / {kT:.5f})) / {Z1:.5f}")
        print(f"Result: n_{n} = {moles_per_level[i]:.5f} moles\n")

    # Final result in the specified format
    final_answer = tuple(moles_per_level)
    print("The ordered set (E_1, E_2, E_3, E_4, E_5) representing the number of moles in each energy level is:")
    print(final_answer)
    print(f"<<<{final_answer}>>>")


solve_particle_distribution()