import numpy as np

def solve_statistical_mechanics():
    """
    Calculates the partition function and the most probable number of particles 
    in each energy level for a given system.
    """
    # Step 1: Define constants
    epsilon_0_meV = 6.9  # meV
    T_K = 4200.0         # Kelvin
    N_total_moles = 3.0  # Total number of particles in moles

    # Physical constants
    k_B_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K
    epsilon_0_eV = epsilon_0_meV * 1e-3 # Convert meV to eV

    # Step 2: Calculate energy levels and degeneracies
    n_levels = np.arange(1, 6, dtype=object) # n = 1, 2, 3, 4, 5. Use dtype=object for large numbers in n^n

    # Calculate epsilon_n for each n
    epsilon_n = epsilon_0_eV * (1 + n_levels**n_levels)

    # Calculate energy levels E_i as cumulative sums of epsilon_n
    E_i = np.cumsum(epsilon_n)

    # Calculate degeneracies g_i
    g_i = 2 * n_levels + 1

    # Step 3: Compute the single-particle partition function (z)
    # The total partition function for N distinguishable particles is Z = z^N.
    # We calculate z, the single-particle partition function.
    
    # Calculate the argument of the exponential, E_i / (k_B * T)
    beta_E = E_i / (k_B_eV_K * T_K)

    # Calculate individual terms for the partition function sum
    z_terms = g_i * np.exp(-beta_E)

    # Sum the terms to get the single-particle partition function z
    z = np.sum(z_terms)

    # Step 4: Calculate the most probable number of particles (in moles) for each level
    # Probability of a particle being in level i
    P_i = z_terms / z

    # Number of moles in each level i
    n_moles_i = N_total_moles * P_i
    
    # Step 5: Output the results
    print(f"The single-particle partition function z is: {z}")
    print("The most probable number of particles in each energy level (in moles) is:")
    
    # Print the final equation with the calculated values
    result_tuple = tuple(n_moles_i)
    print(f"(E1, E2, E3, E4, E5) = ({result_tuple[0]}, {result_tuple[1]}, {result_tuple[2]}, {result_tuple[3]}, {result_tuple[4]})")
    
    # Output the final answer in the specified format
    print(f"<<<{result_tuple}>>>")

solve_statistical_mechanics()