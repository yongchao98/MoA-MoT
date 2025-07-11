import math

def solve_statistical_mechanics():
    """
    Solves the statistical mechanics problem to find the partition function
    and the number of moles in each energy level.
    """
    # 1. Define Constants
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T = 4200.0  # Kelvin
    # Boltzmann constant in eV/K
    k_B = 8.617333262145e-5
    # Total number of moles
    total_moles = 3.0

    # 2. Calculate Energy Levels (En) and Degeneracies (gn)
    num_levels = 5
    epsilons = []
    energies = []
    degeneracies = []
    
    current_energy_sum = 0.0
    for n in range(1, num_levels + 1):
        epsilon_n = epsilon_0_eV * (1 + n**n)
        current_energy_sum += epsilon_n
        
        epsilons.append(epsilon_n)
        energies.append(current_energy_sum)
        
        g_n = 2 * n + 1
        degeneracies.append(g_n)

    # 3. Calculate the Single-Particle Partition Function (Z)
    partition_function = 0.0
    beta = 1.0 / (k_B * T)
    terms = []

    for i in range(num_levels):
        term = degeneracies[i] * math.exp(-energies[i] * beta)
        terms.append(term)
        partition_function += term

    # Print the equation for the partition function with numerical values
    print("The single-particle partition function Z is given by the sum:")
    print("Z = g1*exp(-E1/(k*T)) + g2*exp(-E2/(k*T)) + ... + g5*exp(-E5/(k*T))\n")
    print("Plugging in the values:")
    
    kt_val = k_B * T
    z_eq_str = "Z = "
    for i in range(num_levels):
        g = degeneracies[i]
        E = energies[i]
        z_eq_str += f"{g}*exp(-{E:.4f}/{kt_val:.4f})"
        if i < num_levels - 1:
            z_eq_str += " + "
    print(z_eq_str)

    z_val_str = "Z = "
    for i in range(num_levels):
        z_val_str += f"{terms[i]:.4f}"
        if i < num_levels - 1:
            z_val_str += " + "
    print(z_val_str)
    
    print(f"\nCalculated Partition Function Z = {partition_function:.4f}\n")
    
    # 4. Calculate Occupation Numbers (in Moles)
    moles_in_levels = []
    print("Most probable number of moles in each energy level:")
    for i in range(num_levels):
        # Probability P_n = term_n / Z
        probability_n = terms[i] / partition_function
        # Moles_n = total_moles * P_n
        moles_n = total_moles * probability_n
        moles_in_levels.append(moles_n)
        print(f"Moles in E{i+1}: {moles_n:.4f} moles")

    # 5. Final Output
    final_answer = tuple(moles_in_levels)
    print(f"\nThe ordered set of moles is (E1, E2, E3, E4, E5).")
    
    # Using a f-string to format the final tuple with 4 decimal places for readability
    formatted_answer = f"({', '.join(f'{n:.4f}' for n in final_answer)})"
    print(f"<<<{formatted_answer}>>>")

if __name__ == '__main__':
    solve_statistical_mechanics()