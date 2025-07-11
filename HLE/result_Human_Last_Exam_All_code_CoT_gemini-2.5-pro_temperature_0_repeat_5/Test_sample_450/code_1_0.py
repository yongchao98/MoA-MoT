import math

def solve_statistical_mechanics():
    """
    Calculates the single-particle partition function and the most probable
    number of moles of particles in each of five energy levels.
    """
    # Step 1: Define Constants
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert to eV
    T = 4200.0  # Kelvin
    k_B = 8.617333262e-5  # Boltzmann constant in eV/K
    total_moles = 3.0

    # Calculate thermal energy
    kBT = k_B * T

    # Step 2: Calculate Energy Levels and Degeneracies
    num_levels = 5
    epsilons = []
    energies = []
    degeneracies = []
    
    current_energy = 0.0
    for n in range(1, num_levels + 1):
        # Epsilon_n
        epsilon_n = epsilon_0_eV * (1 + n**n)
        epsilons.append(epsilon_n)
        
        # Cumulative Energy E_n
        current_energy += epsilon_n
        energies.append(current_energy)
        
        # Degeneracy g_n
        g_n = 2 * n + 1
        degeneracies.append(g_n)

    # Step 3: Calculate the Single-Particle Partition Function (Z1)
    Z1 = 0.0
    boltzmann_factors = []
    for i in range(num_levels):
        g_n = degeneracies[i]
        E_n = energies[i]
        # The term in the sum for Z1
        term = g_n * math.exp(-E_n / kBT)
        boltzmann_factors.append(term)
        Z1 += term

    # Step 4: Calculate the Number of Moles in Each Level
    moles_in_levels = []
    for i in range(num_levels):
        # Probability P_n = term / Z1
        probability_n = boltzmann_factors[i] / Z1
        # Moles in level n = total_moles * P_n
        moles_n = total_moles * probability_n
        moles_in_levels.append(moles_n)

    # Step 5: Output the Results
    print(f"The single-particle partition function Z1 is: {Z1}")
    print("\nThe most probable number of particles in each energy level (in moles) is:")
    
    # Create a formatted string for the final tuple
    # The problem asks for the final equation, so we show the sum
    moles_str_list = [f"{mole:.6f}" for mole in moles_in_levels]
    total_moles_check = sum(moles_in_levels)
    
    print(f"({', '.join(moles_str_list)})")
    print(f"\nVerification: The sum of moles is {moles_str_list[0]} + {moles_str_list[1]} + {moles_str_list[2]} + {moles_str_list[3]} + {moles_str_list[4]} = {total_moles_check:.6f}, which equals the total of {total_moles} moles.")
    
    # Final answer in the required format
    final_answer_tuple = tuple(moles_in_levels)
    print(f"\n<<<{final_answer_tuple}>>>")

# Execute the function
solve_statistical_mechanics()