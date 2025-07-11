import numpy as np

def calculate_molar_occupancy():
    """
    Calculates the most probable number of particles (in moles) for each energy level
    in a system of distinguishable particles.
    """
    # 1. Define Constants
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T_K = 4200.0  # Kelvin
    k_B_eV_K = 8.617333262145e-5  # Boltzmann constant in eV/K
    N_total_moles = 3.0 # Total number of moles

    # 2. Calculate Energy Levels and Degeneracies
    num_levels = 5
    energies_eV = []
    degeneracies = []
    
    cumulative_energy_eV = 0.0
    for n in range(1, num_levels + 1):
        # Degeneracy for level n
        g_n = 2 * n + 1
        degeneracies.append(g_n)
        
        # Epsilon for level n
        epsilon_n = epsilon_0_eV * (1 + n**n)
        
        # Cumulative energy E_n
        cumulative_energy_eV += epsilon_n
        energies_eV.append(cumulative_energy_eV)

    # Use numpy arrays for vectorized calculations
    energies_eV = np.array(energies_eV)
    degeneracies = np.array(degeneracies)

    # 3. Calculate the Single-Particle Partition Function (Z1)
    beta = 1.0 / (k_B_eV_K * T_K)
    
    # Calculate terms for the partition function sum: g_j * exp(-beta * E_j)
    z_terms = degeneracies * np.exp(-beta * energies_eV)
    
    # Sum the terms to get Z1
    Z1 = np.sum(z_terms)

    # 4. Calculate the Most Probable Number of Particles (in Moles)
    # Probability P_j = (g_j * exp(-beta * E_j)) / Z1
    # Moles_j = N_total_moles * P_j
    moles_in_level = N_total_moles * z_terms / Z1

    # 5. Output the Result
    # Create the final tuple for printing
    result_tuple = tuple(moles_in_level)
    
    # Print the equation representing the calculation for Z1
    print("Single-particle partition function Z₁ is the sum over levels n=1 to 5:")
    z1_eq_str = "Z₁ = " + " + ".join([f"{g}*exp(-{e:.4f}/(k_B*T))" for g, e in zip(degeneracies, energies_eV)])
    print(z1_eq_str)
    print(f"k_B*T = {k_B_eV_K * T_K:.4f} eV")
    print(f"Z₁ = {Z1:.4f}\n")
    
    print("Most probable number of particles in each level (in moles):")
    # Print the equation for each mole calculation
    for i in range(num_levels):
        # Since Moles_j = 3 * P_j = 3 * z_terms[i] / Z1
        # and z_terms[i] = g_i * exp(-beta * E_i)
        print(f"Moles in E{i+1}: 3 * ({degeneracies[i]} * exp(-{energies_eV[i]:.4f} / {k_B_eV_K * T_K:.4f})) / {Z1:.4f} = {moles_in_level[i]:.4f}")

    print("\nFinal ordered set (E₁, E₂, E₃, E₄, E₅):")
    # The final answer format is just the tuple of numbers
    print(result_tuple)
    
    # Also print the answer in the requested format
    print(f"\n<<<{result_tuple}>>>")

calculate_molar_occupancy()