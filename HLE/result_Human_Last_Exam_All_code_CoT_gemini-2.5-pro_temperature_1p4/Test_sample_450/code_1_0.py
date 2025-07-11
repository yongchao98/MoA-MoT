import numpy as np

def solve_statistical_mechanics():
    """
    Solves for the partition function and the number of moles in each energy level.
    """
    # 1. Define Physical Constants and Parameters
    e0_meV = 6.9  # meV
    e0_eV = e0_meV * 1e-3  # Convert to eV
    T_K = 4200.0  # Kelvin
    k_B_eV_K = 8.617333262e-5  # Boltzmann constant in eV/K
    total_moles = 3.0
    
    # Energy level indices
    n_values = np.arange(1, 6)

    # 2. Calculate Energy Levels (E_n) and Degeneracies (g_n)
    g_n = 2 * n_values + 1
    eps_n_eV = e0_eV * (1 + n_values**n_values)
    E_n_eV = np.cumsum(eps_n_eV)
    
    k_B_T_eV = k_B_eV_K * T_K

    print("--- Calculating Energy Levels and Degeneracies ---")
    for i in range(len(n_values)):
        print(f"Level n={n_values[i]}: E_{n_values[i]} = {E_n_eV[i]:.4f} eV, g_{n_values[i]} = {g_n[i]}")
    print(f"Thermal energy k_B*T = {k_B_T_eV:.4f} eV\n")
    
    # 3. Calculate the Single-Particle Partition Function (Z_1)
    print("--- Calculating the Single-Particle Partition Function (Z_1) ---")
    print("Z_1 = sum(g_n * exp(-E_n / (k_B * T)))")

    # Constructing the equation string for Z_1
    z1_eq_str = "Z_1 = " + " + ".join([f"{g}*exp(-{E:.4f}/{k_B_T_eV:.4f})" for g, E in zip(g_n, E_n_eV)])
    print(z1_eq_str)
    
    boltzmann_terms = g_n * np.exp(-E_n_eV / k_B_T_eV)
    Z_1 = np.sum(boltzmann_terms)

    z1_val_str = "Z_1 = " + " + ".join([f"{term:.4f}" for term in boltzmann_terms]) + f" = {Z_1:.4f}"
    print(z1_val_str)
    print("Note: The total partition function Z = (Z_1)^N, where N is the total number of particles. Z_1 is the single-particle partition function.\n")

    # 4. Calculate the Most Probable Number of Moles in Each Energy Level (n_n)
    print("--- Calculating Number of Moles in Each Energy Level ---")
    
    mole_numbers = []
    for i in range(len(n_values)):
        n_val = n_values[i]
        g_val = g_n[i]
        E_val = E_n_eV[i]
        boltzmann_term = boltzmann_terms[i]

        # Calculate moles for the current level
        moles = total_moles * boltzmann_term / Z_1
        mole_numbers.append(moles)
        
        print(f"For level n={n_val}:")
        print(f"  n_{n_val} = N_moles * (g_{n_val} * exp(-E_{n_val} / (k_B*T))) / Z_1")
        print(f"  n_{n_val} = {total_moles} * ({g_val} * exp(-{E_val:.4f} / {k_B_T_eV:.4f})) / {Z_1:.4f}")
        print(f"  n_{n_val} = {total_moles} * {boltzmann_term:.4f} / {Z_1:.4f} = {moles:.4f} moles")

    # 5. Output the final result in the specified format
    final_answer = tuple(mole_numbers)
    print("\nThe ordered set representing the number of moles in each energy level (E_1, E_2, E_3, E_4, E_5) is:")
    print(f"{final_answer}")
    
    # Final answer in the required format
    print(f"\n<<<{final_answer}>>>")


solve_statistical_mechanics()