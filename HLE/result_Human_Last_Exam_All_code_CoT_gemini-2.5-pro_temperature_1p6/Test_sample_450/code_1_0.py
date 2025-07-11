import math

def solve_statistical_mechanics():
    """
    Solves the statistical mechanics problem to find the partition function
    and the most probable number of moles in each energy level.
    """
    # 1. Define constants and parameters
    epsilon_0_meV = 6.9
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    k_B_eV_K = 8.617333262e-5           # Boltzmann constant in eV/K
    T_K = 4200.0                         # Temperature in Kelvin
    N_total_moles = 3.0                  # Total number of moles of particles
    num_levels = 5

    # Pre-calculate the thermal energy
    k_B_T_eV = k_B_eV_K * T_K

    # 2. Calculate energy levels and degeneracies
    energies = []
    degeneracies = []
    cumulative_energy = 0.0

    for n in range(1, num_levels + 1):
        # Calculate epsilon_n
        epsilon_n = epsilon_0_eV * (1 + n**n)
        
        # Calculate cumulative energy E_n
        cumulative_energy += epsilon_n
        energies.append(cumulative_energy)

        # Calculate degeneracy g_n
        g_n = 2 * n + 1
        degeneracies.append(g_n)

    # 3. Calculate the single-particle partition function (Z)
    partition_Z = 0.0
    boltzmann_terms = []
    for i in range(num_levels):
        term = degeneracies[i] * math.exp(-energies[i] / k_B_T_eV)
        boltzmann_terms.append(term)
        partition_Z += term

    print(f"Calculated single-particle partition function Z = {partition_Z:.5f}")
    print(f"Thermal energy k_B*T = {k_B_T_eV:.5f} eV")
    print("-" * 50)
    print("The final equation for the number of moles in each level is:")
    print("moles_n = (Total Moles) * (g_n * exp(-E_n / k_B*T)) / Z")
    print("-" * 50)


    # 4. Calculate and display the number of moles in each level
    moles_per_level = []
    for i in range(num_levels):
        g_n = degeneracies[i]
        E_n = energies[i]
        
        # Calculate moles for the current level
        moles_n = N_total_moles * boltzmann_terms[i] / partition_Z
        moles_per_level.append(moles_n)
        
        # Print the detailed equation with values for this level
        print(f"moles_{i+1} = {N_total_moles} * ({g_n} * exp(-{E_n:.5f} eV / {k_B_T_eV:.5f} eV)) / {partition_Z:.5f} = {moles_n:.5f} moles")

    # Final sanity check
    print("-" * 50)
    print(f"Check: Total moles = {sum(moles_per_level):.5f} (should be close to {N_total_moles})")
    
    # Format final answer
    final_answer_tuple = tuple(moles_per_level)
    
    # Using a f-string to format the tuple with 5 decimal places for clarity
    formatted_tuple_string = '(' + ', '.join(f'{x:.5f}' for x in final_answer_tuple) + ')'
    
    print(f"\nFinal Answer: The ordered set of moles in each energy level is {formatted_tuple_string}")

    # Return the answer in the required format
    print(f"\n<<<{formatted_tuple_string}>>>")


solve_statistical_mechanics()