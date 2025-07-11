import numpy as np

def solve_statistical_mechanics_problem():
    """
    Calculates the single-particle partition function and the number of moles
    in each of five energy levels for a system of distinguishable particles.
    """
    # Step 1: Define Constants
    epsilon_0_meV = 6.9
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T_K = 4200.0  # Temperature in Kelvin
    k_B_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K
    total_moles = 3.0

    # Thermal energy
    kT_eV = k_B_eV_K * T_K

    # Step 2: Calculate Energy Levels and Degeneracies
    num_levels = 5
    epsilons = []
    energies = []
    degeneracies = []
    
    current_energy = 0.0
    for n in range(1, num_levels + 1):
        # Calculate epsilon_n for the current level
        epsilon_n = epsilon_0_eV * (1 + n**n)
        epsilons.append(epsilon_n)
        
        # Energy levels are cumulative sums
        current_energy += epsilon_n
        energies.append(current_energy)
        
        # Calculate degeneracy
        g_n = 2 * n + 1
        degeneracies.append(g_n)

    # Step 3: Calculate the Single-Particle Partition Function (Z1)
    Z1 = 0.0
    boltzmann_terms = []
    for i in range(num_levels):
        g_n = degeneracies[i]
        E_n = energies[i]
        exponent = -E_n / kT_eV
        term = g_n * np.exp(exponent)
        boltzmann_terms.append(term)
        Z1 += term

    # Step 4: Calculate Moles per Energy Level
    moles_per_level = []
    for term in boltzmann_terms:
        # Probability of a particle being in the level
        P_n = term / Z1
        # Number of moles in the level
        n_moles_n = total_moles * P_n
        moles_per_level.append(n_moles_n)

    # Step 5: Output the Results
    print(f"Physical Constants:")
    print(f"epsilon_0 = {epsilon_0_eV:.5f} eV")
    print(f"T = {T_K} K")
    print(f"k_B*T = {kT_eV:.5f} eV\n")
    
    print("Calculated Energy Levels and Degeneracies:")
    for i in range(num_levels):
        print(f"Level n={i+1}: E_{i+1} = {energies[i]:.5f} eV, g_{i+1} = {degeneracies[i]}")

    print(f"\nSingle-Particle Partition Function Z_1 = {Z1:.5f}\n")

    print("Most probable number of particles in each energy level (in moles):")
    # The problem asks to output each number in the final equation.
    # The final equation is the tuple of moles.
    result_tuple = tuple(moles_per_level)
    print(f"({result_tuple[0]:.5f}, {result_tuple[1]:.5f}, {result_tuple[2]:.5f}, {result_tuple[3]:.5f}, {result_tuple[4]:.5f})")
    
    # Final answer in the required format
    # Formatting to 5 decimal places for consistency
    final_answer_str = f"({result_tuple[0]:.5f}, {result_tuple[1]:.5f}, {result_tuple[2]:.5f}, {result_tuple[3]:.5f}, {result_tuple[4]:.5f})"
    print(f"\n<<<{final_answer_str}>>>")


if __name__ == "__main__":
    solve_statistical_mechanics_problem()