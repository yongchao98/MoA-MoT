import math

def solve_statistical_mechanics():
    """
    Calculates the partition function and the most probable number of particles
    in each energy level for a given system.
    """
    # Step 1: Define Constants
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T = 4200.0  # Kelvin
    k_B_eV_K = 8.617333262145e-5  # Boltzmann constant in eV/K
    total_moles = 3.0
    num_levels = 5

    # Thermal energy
    kBT = k_B_eV_K * T

    # Lists to store calculated values
    epsilons = []
    energies_E_n = []
    degeneracies_g_n = []

    # Step 2: Calculate Energy Levels (E_n)
    # First, calculate epsilon_n increments
    for n in range(1, num_levels + 1):
        epsilon_n = epsilon_0_eV * (1 + n**n)
        epsilons.append(epsilon_n)
    
    # Then, calculate the cumulative energy E_n
    current_energy = 0.0
    for epsilon in epsilons:
        current_energy += epsilon
        energies_E_n.append(current_energy)

    # Step 3: Calculate Degeneracies (g_n)
    for n in range(1, num_levels + 1):
        g_n = 2 * n + 1
        degeneracies_g_n.append(g_n)

    # Step 4: Calculate Single-Particle Partition Function (Z1)
    Z1_terms = []
    for i in range(num_levels):
        E_n = energies_E_n[i]
        g_n = degeneracies_g_n[i]
        # Boltzmann factor term for the partition function sum
        term = g_n * math.exp(-E_n / kBT)
        Z1_terms.append(term)
    
    Z1 = sum(Z1_terms)

    # Step 5: Calculate Population in Moles
    moles_per_level = []
    for term in Z1_terms:
        # Probability P_n = term / Z1
        # Moles_n = total_moles * P_n
        moles = total_moles * (term / Z1)
        moles_per_level.append(moles)

    # Step 6: Format the Output
    print(f"The thermal energy k_B * T is: {kBT:.4f} eV")
    print(f"The single-particle partition function Z1 is: {Z1:.4f}")
    print("\n--- Most Probable Number of Moles per Energy Level ---")
    
    # The final equation for moles is: moles_n = total_moles * (g_n * exp(-E_n / kBT)) / Z1
    # We print each component for clarity
    for i in range(num_levels):
        print(f"\nLevel E{i+1}:")
        print(f"  Energy E{i+1} = {energies_E_n[i]:.4f} eV")
        print(f"  Degeneracy g{i+1} = {degeneracies_g_n[i]}")
        print(f"  Boltzmann term (g{i+1}*exp(-E{i+1}/kBT)) = {Z1_terms[i]:.4f}")
        final_moles = moles_per_level[i]
        print(f"  Final Equation: {total_moles:.1f} * ({Z1_terms[i]:.4f} / {Z1:.4f}) = {final_moles:.4f} moles")

    # Final answer in the specified format
    final_answer_tuple = tuple(moles_per_level)
    print("\nFinal ordered set for moles in levels (E1, E2, E3, E4, E5):")
    print(f"<<<{final_answer_tuple}>>>")

# Execute the function
solve_statistical_mechanics()