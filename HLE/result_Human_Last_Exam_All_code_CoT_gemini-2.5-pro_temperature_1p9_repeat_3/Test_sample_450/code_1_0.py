import math

def calculate_particle_distribution():
    """
    Calculates the partition function and the most probable number of moles
    in each energy level for a system of distinguishable particles.
    """
    # 1. Define constants and parameters
    e0_meV = 6.9
    T_K = 4200.0
    total_moles = 3.0
    
    # Convert energy to eV
    e0_eV = e0_meV * 1e-3
    
    # Physical constant: Boltzmann constant in eV/K
    kB_eV_per_K = 8.617333262e-5
    
    # Thermal energy in eV
    kT = kB_eV_per_K * T_K
    
    # 2. Calculate energy levels and degeneracies
    num_levels = 5
    energies_eV = []
    degeneracies = []
    
    current_energy_sum = 0.0
    for n in range(1, num_levels + 1):
        # Epsilon component: epsilon_n(n) = e0 * (1 + n^n)
        epsilon_n = e0_eV * (1 + n**n)
        
        # Energy level E_n is the cumulative sum of epsilons
        current_energy_sum += epsilon_n
        energies_eV.append(current_energy_sum)
        
        # Degeneracy: g_n = 2n + 1
        g_n = 2 * n + 1
        degeneracies.append(g_n)
        
    # 3. Compute the single-particle partition function Z
    Z = 0.0
    boltzmann_terms = []
    for i in range(num_levels):
        E_n = energies_eV[i]
        g_n = degeneracies[i]
        
        # Calculate the term g_n * exp(-E_n / kT) for the sum
        term = g_n * math.exp(-E_n / kT)
        boltzmann_terms.append(term)
        Z += term
        
    # 4. Calculate the most probable number of moles in each level
    moles_in_levels = []
    for term in boltzmann_terms:
        # Probability P_n = (term / Z)
        # Moles n_n = total_moles * P_n
        moles_n = total_moles * (term / Z)
        moles_in_levels.append(moles_n)
        
    # 5. Print the results
    print(f"The single-particle partition function Z is: {Z:.4f}")
    
    # The problem asks for the ordered set representing moles in levels E1, E2, ...
    print("\nThe most probable number of moles in each energy level (n_E1, n_E2, n_E3, n_E4, n_E5) is:")
    
    # Constructing the final output string showing each number in the set
    result_string = f"({moles_in_levels[0]:.4f}, {moles_in_levels[1]:.4f}, {moles_in_levels[2]:.4f}, {moles_in_levels[3]:.4f}, {moles_in_levels[4]:.4f})"
    print(result_string)

# Execute the function
calculate_particle_distribution()
