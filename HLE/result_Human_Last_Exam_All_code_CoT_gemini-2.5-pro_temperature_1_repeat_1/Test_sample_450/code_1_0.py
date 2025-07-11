import math

def calculate_particle_distribution():
    """
    Calculates the most probable number of particles (in moles) in each of five energy levels
    for a system of distinguishable particles in thermal equilibrium.
    """
    # 1. Define Constants and Parameters
    e0_meV = 6.9  # meV
    e0_eV = e0_meV * 1e-3  # Convert meV to eV
    T_K = 4200  # Kelvin
    N_moles_total = 3.0  # Total number of moles
    kB_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K

    # Pre-calculate thermal energy
    kBT = kB_eV_K * T_K

    # Lists to store values for each level
    energies = []
    degeneracies = []
    boltzmann_terms = []
    
    # 2. Calculate Energy Levels and Degeneracies
    current_E = 0.0
    print("Calculating properties for each energy level (n=1 to 5):")
    for n in range(1, 6):
        # Calculate energy contribution for level n
        epsilon_n = e0_eV * (1 + n**n)
        
        # Energy levels are cumulative
        current_E += epsilon_n
        energies.append(current_E)
        
        # Calculate degeneracy
        g_n = 2 * n + 1
        degeneracies.append(g_n)
        
        # Calculate the term for the partition function sum
        term = g_n * math.exp(-current_E / kBT)
        boltzmann_terms.append(term)
        
        print(f"  Level n={n}: E_{n} = {current_E:.4f} eV, g_{n} = {g_n}")

    # 3. Compute the Single-Particle Partition Function (Z1)
    Z1 = sum(boltzmann_terms)
    print(f"\nThe single-particle partition function Z1 is: {Z1:.4f}")

    # 4. Calculate the Number of Moles in Each Level
    moles_in_levels = []
    for term in boltzmann_terms:
        # Probability Pn = term / Z1
        # Moles = N_total_moles * Pn
        moles_n = N_moles_total * (term / Z1)
        moles_in_levels.append(moles_n)

    # 5. Output the Result
    print("\nThe most probable number of particles in each energy level, in moles, is:")
    print("(moles in E1, moles in E2, moles in E3, moles in E4, moles in E5)")
    
    # Format the final tuple for clean output
    result_tuple = tuple(moles_in_levels)
    print(f"Final Answer: {result_tuple}")

# Run the calculation
calculate_particle_distribution()

<<< (0.7963690623259441, 1.2064104439058564, 0.9880196238699198, 0.00919936993181881, 1.5000033539121648e-26) >>>