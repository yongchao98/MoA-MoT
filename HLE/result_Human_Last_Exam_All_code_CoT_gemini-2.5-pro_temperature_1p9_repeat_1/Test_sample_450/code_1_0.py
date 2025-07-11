import math

def solve_particle_distribution():
    """
    Calculates the most probable number of particles in each energy level for a given system.
    """
    # 1. Define Constants
    epsilon_0_meV = 6.9
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T_K = 4200.0  # Temperature in Kelvin
    k_B_eV_K = 8.617333e-5  # Boltzmann constant in eV/K
    N_total_moles = 3.0 # Total number of particles in moles

    # Calculate thermal energy
    kbT = k_B_eV_K * T_K

    # Lists to store values for the 5 levels
    energy_levels = [0] * 5
    degeneracies = [0] * 5
    epsilons = [0] * 5

    # 2. Calculate incremental energies (eps_n) and degeneracies (g_n)
    for n in range(1, 6):
        # Degeneracy
        g_n = 2 * n + 1
        degeneracies[n-1] = g_n
        
        # Incremental energy
        eps_n = epsilon_0_eV * (1 + n**n)
        epsilons[n-1] = eps_n

    # 3. Calculate cumulative energy levels (E_n)
    cumulative_energy = 0
    for i in range(5):
        cumulative_energy += epsilons[i]
        energy_levels[i] = cumulative_energy

    # 4. Calculate the single-particle partition function (z)
    partition_function_z = 0
    boltzmann_terms = [0] * 5
    print("Calculating the single-particle partition function z...")
    print("-" * 50)
    print(f"Thermal energy k_B*T = {kbT:.4f} eV")
    for i in range(5):
        n = i + 1
        g_n = degeneracies[i]
        E_n = energy_levels[i]
        boltzmann_factor = math.exp(-E_n / kbT)
        term = g_n * boltzmann_factor
        boltzmann_terms[i] = term
        partition_function_z += term
        print(f"Level n={n}: g_{n}={g_n}, E_{n}={E_n:.4f} eV, Term = {g_n}*exp(-{E_n:.4f}/{kbT:.4f}) = {term:.4f}")

    print("-" * 50)
    print(f"Total single-particle partition function z = {partition_function_z:.4f}")
    print("\n" + "="*50 + "\n")
    
    # 5. Calculate the number of moles in each energy level
    moles_in_levels = [0] * 5
    print("Calculating the number of moles in each energy level (n_i)...")
    print("-" * 50)

    for i in range(5):
        n = i + 1
        term = boltzmann_terms[i]
        moles = N_total_moles * term / partition_function_z
        moles_in_levels[i] = moles
        
        # Output each number in the final equation
        print(f"n_{n} = N_total * (Term_{n} / z)")
        print(f"n_{n} = {N_total_moles} * ({term:.4f} / {partition_function_z:.4f}) = {moles:.4f} moles")
        print("-" * 20)
    
    print("\n" + "="*50 + "\n")

    # 6. Present the final result
    final_moles_tuple = tuple(moles_in_levels)
    print("The most probable number of particles in each energy level, expressed in moles, is:")
    print(f"(n1, n2, n3, n4, n5) = {final_moles_tuple}")
    
    # Final answer in the required format
    answer_string = f"({moles_in_levels[0]:.4f}, {moles_in_levels[1]:.4f}, {moles_in_levels[2]:.4f}, {moles_in_levels[3]:.4f}, {moles_in_levels[4]:.4f})"
    print(f"\n<<<{answer_string}>>>")


solve_particle_distribution()