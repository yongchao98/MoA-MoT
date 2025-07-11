import math

def solve_statistical_mechanics():
    """
    Calculates the partition function and the most probable number of moles
    of particles in each of five energy levels.
    """
    # Step 1: Define Constants
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV / 1000.0  # eV
    e_charge = 1.60218e-19  # C, elementary charge
    epsilon_0 = epsilon_0_eV * e_charge  # J

    T = 4200.0  # K
    k_B = 1.380649e-23  # J/K
    total_moles = 3.0

    # Thermal energy
    k_B_T = k_B * T

    print("--- System Parameters ---")
    print(f"epsilon_0 = {epsilon_0_meV} meV")
    print(f"Temperature T = {T} K")
    print(f"Thermal Energy k_B*T = {k_B_T / e_charge:.4f} eV\n")

    # Step 2: Calculate Energy Levels and Degeneracies
    n_values = range(1, 6)
    epsilons = [epsilon_0 * (1 + n**n) for n in n_values]
    
    energies = []
    current_energy = 0
    for eps in epsilons:
        current_energy += eps
        energies.append(current_energy)
        
    degeneracies = [2 * n + 1 for n in n_values]

    print("--- Energy Levels and Degeneracies ---")
    print("Level (n) | Energy E_n (meV) | Degeneracy g_n")
    print("---------------------------------------------")
    for i, n in enumerate(n_values):
        print(f"{n:^9} | {energies[i] / e_charge * 1000:^18.2f} | {degeneracies[i]:^14}")
    print("\n")

    # Step 3: Calculate the Single-Particle Partition Function (z)
    z_terms = []
    for i in range(len(n_values)):
        exponent = -energies[i] / k_B_T
        term = degeneracies[i] * math.exp(exponent)
        z_terms.append(term)
    
    z = sum(z_terms)

    print("--- Partition Function Calculation ---")
    print("The single-particle partition function z is the sum of terms g_n * exp(-E_n / (k_B*T))")
    for i, n in enumerate(n_values):
        print(f"Term for n={n}: {z_terms[i]:.4f}")
    print(f"\nTotal Partition Function z = {z:.4f}\n")
    print("(Note: The total partition function for N distinguishable particles is Z = z^N)\n")

    # Step 4: Calculate the Number of Moles per Energy Level (n_n)
    moles_per_level = []
    print("--- Moles per Energy Level Calculation ---")
    print("The number of moles in level n, n_n, is given by: 3 * (g_n * exp(-E_n / k_B*T)) / z")

    for i, n in enumerate(n_values):
        moles_n = total_moles * z_terms[i] / z
        moles_per_level.append(moles_n)
        print(f"\nFor n = {n}:")
        print(f"n_{n} = {total_moles:.1f} moles * ({degeneracies[i]} * exp(-{energies[i]/e_charge:.4f} eV / {k_B_T/e_charge:.4f} eV)) / {z:.4f}")
        print(f"n_{n} = {total_moles:.1f} * ({z_terms[i]:.4f}) / {z:.4f}")
        print(f"n_{n} = {moles_n:.4f} moles")

    # Step 5: Output the Final Result
    final_result = tuple(moles_per_level)
    print("\n--- Final Result ---")
    print("The ordered set (E1, E2, E3, E4, E5) representing the number of moles in each energy level is:")
    print(f"({final_result[0]:.4f}, {final_result[1]:.4f}, {final_result[2]:.4f}, {final_result[3]:.4f}, {final_result[4]:.4f})")
    
    # Required final answer format
    final_answer_string = f"<<<({final_result[0]:.4f}, {final_result[1]:.4f}, {final_result[2]:.4f}, {final_result[3]:.4f}, {final_result[4]:.4f})>>>"
    print(final_answer_string)


solve_statistical_mechanics()