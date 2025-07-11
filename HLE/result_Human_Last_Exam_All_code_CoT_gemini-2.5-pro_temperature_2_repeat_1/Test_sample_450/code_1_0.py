import math

def solve_particle_distribution():
    """
    Calculates the single-particle partition function and the most probable
    number of particles in moles for each energy level.
    """
    # Step 1: Define Constants and Parameters
    epsilon_0_meV = 6.9
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T_K = 4200.0  # Temperature in Kelvin
    k_B_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K
    N_moles = 3.0 # Total number of moles of particles

    # Thermal energy
    k_B_T = k_B_eV_K * T_K

    # Lists to store values for each level n=1 to 5
    energies = []
    degeneracies = []
    
    # Step 2: Calculate Energy Levels and Degeneracies
    current_E = 0.0
    boltzmann_terms = []

    print("Calculating energy levels, degeneracies, and Boltzmann terms:\n")
    for n in range(1, 6):
        # Degeneracy
        g_n = 2 * n + 1
        degeneracies.append(g_n)

        # Energy increment
        epsilon_n = epsilon_0_eV * (1 + n**n)
        
        # Cumulative energy level
        current_E += epsilon_n
        energies.append(current_E)
        
        # Boltzmann factor for the partition function sum
        boltzmann_term = g_n * math.exp(-current_E / k_B_T)
        boltzmann_terms.append(boltzmann_term)
        
        print(f"Level n={n}:")
        print(f"  Degeneracy g_{n} = {g_n}")
        print(f"  Energy E_{n} = {current_E:.4f} eV")
        print(f"  Boltzmann Term (g_{n}*exp(-E_{n}/kBT)) = {boltzmann_term:.4f}")
        print("-" * 20)

    # Step 3: Compute the Single-Particle Partition Function (Z_1)
    Z_1 = sum(boltzmann_terms)
    print(f"\nSingle-particle partition function Z_1 = {Z_1:.4f}\n")

    # Step 4: Determine Molar Populations
    molar_populations = []
    print("Calculating the most probable number of moles in each level:\n")
    for i in range(5):
        n = i + 1
        # Probability P_n = (Boltzmann Term)_n / Z_1
        # Moles n_n = N_moles * P_n
        population_moles = N_moles * boltzmann_terms[i] / Z_1
        molar_populations.append(population_moles)
        print(f"Moles in Level E_{n} = {N_moles:.1f} * ({boltzmann_terms[i]:.4f} / {Z_1:.4f}) = {population_moles:.4f} moles")
    
    # Step 5: Output the Final Result
    # Create the string representation of the final ordered set
    final_answer_str = "(" + ", ".join(f"{p:.4f}" for p in molar_populations) + ")"
    
    print("\nThe ordered set representing the number of moles in each energy level (E1, E2, E3, E4, E5) is:")
    print(final_answer_str)

    # Final answer in the required format
    final_answer_tuple = tuple(molar_populations)
    return final_answer_tuple

# Run the calculation and print the final result in the specified format
final_result = solve_particle_distribution()
final_result_formatted_string = "(" + ", ".join(f"{x:.4f}" for x in final_result) + ")"
print(f"\n<<<{final_result_formatted_string}>>>")
