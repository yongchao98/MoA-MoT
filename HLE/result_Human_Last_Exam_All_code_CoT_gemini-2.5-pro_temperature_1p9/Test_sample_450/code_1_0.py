import math

def solve_statistical_mechanics():
    """
    Solves the statistical mechanics problem to find the partition function
    and the distribution of particles across energy levels.
    """
    # 1. Define Physical Constants
    e0_meV = 6.9  # meV
    e0_eV = e0_meV * 1e-3  # Convert meV to eV
    T_K = 4200.0  # Temperature in Kelvin
    k_B_eV_K = 8.617333e-5  # Boltzmann constant in eV/K
    total_moles = 3.0

    # Calculate thermal energy
    kT_eV = k_B_eV_K * T_K

    # 2. Determine Energy Levels and Degeneracies
    num_levels = 5
    energies_eV = []
    degeneracies = []
    
    current_E_eV = 0.0
    for n in range(1, num_levels + 1):
        # Calculate energy increment epsilon_n
        epsilon_n_eV = e0_eV * (1 + n**n)
        
        # Energy levels are cumulative sums
        current_E_eV += epsilon_n_eV
        energies_eV.append(current_E_eV)
        
        # Calculate degeneracy g_n
        g_n = 2 * n + 1
        degeneracies.append(g_n)

    # 3. Calculate the Single-Particle Partition Function (Z1)
    Z1 = 0.0
    terms = []
    for i in range(num_levels):
        g = degeneracies[i]
        E = energies_eV[i]
        # The exponent is -E / kT
        exponent = -E / kT_eV
        term = g * math.exp(exponent)
        terms.append(term)
        Z1 += term
        
    print(f"The single-particle partition function Z1 is: {Z1:.4f}\n")

    # 4. Calculate and display the number of moles in each energy level
    moles_distribution = []
    print("Calculation of moles in each energy level:")
    for i in range(num_levels):
        g = degeneracies[i]
        E = energies_eV[i]
        
        # Probability for a particle to be in level i
        prob_i = terms[i] / Z1
        
        # Number of moles in level i
        moles_i = total_moles * prob_i
        moles_distribution.append(moles_i)

        # Print the equation with numerical values as requested
        # Moles(E_i) = Total Moles * (g_i * exp(-E_i / kT)) / Z1
        print(f"Moles(E{i+1}) = {total_moles:.1f} mol * ({g} * exp(-{E:.4f} eV / {kT_eV:.4f} eV)) / {Z1:.4f} = {moles_i:.4f} mol")

    # 5. Output the final answer in the specified format
    final_answer = tuple(moles_distribution)
    print(f"\nThe ordered set for the number of moles in each energy level (E1, E2, E3, E4, E5) is:")
    print(f"<<<{final_answer}>>>")

# Execute the function to get the solution
solve_statistical_mechanics()