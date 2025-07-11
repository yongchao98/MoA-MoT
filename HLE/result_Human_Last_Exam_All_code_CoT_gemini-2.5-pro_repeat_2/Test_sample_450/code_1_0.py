import numpy as np

def solve_particle_distribution():
    """
    Calculates the partition function and the most probable number of moles
    of particles in each of five specified energy levels.
    """
    # 1. Define Constants
    epsilon_0_meV = 6.9
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T_K = 4200.0  # Temperature in Kelvin
    k_B_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K
    N_moles_total = 3.0  # Total number of moles

    # Calculate thermal energy for convenience
    kT = k_B_eV_K * T_K

    # 2. Calculate Energy Levels and Degeneracies
    n_levels = 5
    n_values = np.arange(1, n_levels + 1)

    # Degeneracies: g_n = 2n + 1
    g_n = 2 * n_values + 1

    # Epsilon values: epsilon_n = epsilon_0 * (1 + n^n)
    # Note: n**n can lead to very large numbers, ensure float type for calculations
    epsilon_n = epsilon_0_eV * (1 + n_values.astype(float)**n_values.astype(float))

    # Energy levels E_n are the cumulative sum of epsilon_n
    E_n = np.cumsum(epsilon_n)

    # 3. Compute the Partition Function (Z)
    # Boltzmann factors for each level
    boltzmann_factors = np.exp(-E_n / kT)

    # Each term in the partition function sum: g_n * exp(-E_n/kT)
    z_terms = g_n * boltzmann_factors
    
    # The single-particle partition function Z is the sum of these terms
    Z = np.sum(z_terms)

    # 4. Calculate Moles per Level
    # Probability of a particle being in level n
    probabilities = z_terms / Z
    
    # Most probable number of moles in each level
    n_moles_levels = N_moles_total * probabilities

    # 5. Format and Print Results
    print(f"The thermal energy kT is: {kT:.5f} eV")
    print(f"The single-particle partition function Z is: {Z:.5f}")
    print("\n--- Calculation of Moles per Energy Level ---\n")

    final_moles_list = []
    for i in range(n_levels):
        moles = n_moles_levels[i]
        final_moles_list.append(moles)
        
        print(f"Level E{i+1}:")
        print(f"  - Energy E_{i+1} = {E_n[i]:.5f} eV")
        print(f"  - Degeneracy g_{i+1} = {g_n[i]}")
        print(f"  - Boltzmann Term g_{i+1}*exp(-E_{i+1}/kT) = {z_terms[i]:.5f}")
        print(f"  - Probability P_{i+1} = {z_terms[i]:.5f} / {Z:.5f} = {probabilities[i]:.5f}")
        print(f"  - Equation for moles: n_{i+1} = N_total * P_{i+1}")
        print(f"  - Number of Moles n_{i+1} = {N_moles_total:.1f} * {probabilities[i]:.5f} = {moles:.5f} moles\n")

    # Format the final answer as an ordered set
    result_tuple_str = "(" + ", ".join([f"{n:.5f}" for n in final_moles_list]) + ")"
    
    print("--------------------------------------------------")
    print(f"The final ordered set of moles (n1, n2, n3, n4, n5) is:")
    print(result_tuple_str)
    print("--------------------------------------------------")
    
    # Final answer in the required format
    print(f"<<<{result_tuple_str}>>>")

# Execute the function
solve_particle_distribution()