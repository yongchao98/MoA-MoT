import math

def calculate_particle_distribution():
    """
    Calculates the most probable distribution of particles in moles across five energy levels.
    """
    # Step 1: Define Constants
    epsilon_0_meV = 6.9  # meV
    eV_to_J = 1.60218e-19  # J/eV
    epsilon_0 = epsilon_0_meV * 1e-3 * eV_to_J  # Convert meV to Joules
    
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    T = 4200.0  # Temperature in Kelvin
    
    # Total number of moles of particles
    total_moles = 3.0
    
    # Thermal energy
    k_B_T = k_B * T

    # Step 2: Calculate energy levels (E_n) and degeneracies (g_n)
    num_levels = 5
    energies = []
    degeneracies = []
    
    cumulative_energy = 0.0
    for n in range(1, num_levels + 1):
        # Calculate epsilon_n
        epsilon_n = epsilon_0 * (1 + n**n)
        
        # Calculate cumulative energy E_n
        cumulative_energy += epsilon_n
        energies.append(cumulative_energy)
        
        # Calculate degeneracy g_n
        degeneracy = 2 * n + 1
        degeneracies.append(degeneracy)

    # Step 3: Calculate the single-particle partition function (Z)
    partition_function_Z = 0.0
    boltzmann_terms = []
    
    for i in range(num_levels):
        E_n = energies[i]
        g_n = degeneracies[i]
        
        # Calculate the Boltzmann factor term for the partition function
        term = g_n * math.exp(-E_n / k_B_T)
        boltzmann_terms.append(term)
        partition_function_Z += term
        
    # Step 4: Calculate the most probable number of moles in each level
    moles_distribution = []
    for term in boltzmann_terms:
        # Probability of occupying level n
        probability = term / partition_function_Z
        # Number of moles in level n
        moles_n = total_moles * probability
        moles_distribution.append(moles_n)
        
    # Step 5: Print the result
    # The phrase "final equation" is interpreted as showing the final calculated values.
    # We present the result as a tuple, which is an ordered set of numbers.
    print(f"The single-particle partition function Z is: {partition_function_Z}")
    print("The most probable number of particles in moles for each energy level (E1, E2, E3, E4, E5) is:")
    
    # Printing each number in the final tuple
    result_tuple = tuple(moles_distribution)
    # The format (val1, val2, ...) represents the ordered set as requested.
    # Printing each individual number of the final tuple, as per the instruction.
    print(f"({result_tuple[0]}, {result_tuple[1]}, {result_tuple[2]}, {result_tuple[3]}, {result_tuple[4]})")


# Execute the calculation and print the results
calculate_particle_distribution()

# The final answer in the requested format. I am formatting the tuple to 5 decimal places.
# moles_distribution = [0.79549, 1.20524, 0.98971, 0.00953, 2.76e-26]
# result = (0.79549, 1.20524, 0.98971, 0.00953, 0.00000)
# Wrapping this in the requested <<<...>>> format.
<<<(0.79549, 1.20524, 0.98971, 0.00953, 0.00000)>>>