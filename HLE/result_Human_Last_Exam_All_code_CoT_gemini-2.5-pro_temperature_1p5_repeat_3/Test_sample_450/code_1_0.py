import numpy as np

def calculate_particle_distribution():
    """
    Calculates the single-particle partition function and the most probable number of moles
    of particles in each of five specified energy levels.
    """
    # 1. Define Constants
    epsilon_0_meV = 6.9  # meV
    epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
    T = 4200.0  # Kelvin
    k_B = 8.617333262145e-5  # Boltzmann constant in eV/K
    total_moles = 3.0  # Total moles of particles (N = 3 * N_A)

    # Thermal energy
    kT = k_B * T

    # Level indices
    n = np.arange(1, 6)

    # 2. Calculate Energy Levels and Degeneracies
    # Calculate energy shifts epsilon_n(n) = epsilon_0 * (1 + n^n)
    epsilon_n = epsilon_0_eV * (1 + n**n)

    # Calculate energy levels E_n as the cumulative sum of epsilon_n
    E_n = np.cumsum(epsilon_n)

    # Calculate degeneracies g_n = 2n + 1
    g_n = 2 * n + 1

    # 3. Calculate the Single-Particle Partition Function (Z1)
    # The particles are distinguishable, so the relevant quantity for population calculation
    # is the single-particle partition function, Z1.
    boltzmann_factors = g_n * np.exp(-E_n / kT)
    Z1 = np.sum(boltzmann_factors)

    # 4. Calculate Level Populations
    # Probability of a particle being in level n
    P_n = boltzmann_factors / Z1
    
    # Most probable number of moles in each level
    moles_in_level = total_moles * P_n

    # 5. Output Results
    print(f"The single-particle partition function Z1 is: {Z1}")
    
    print("\nThe most probable number of particles in each energy level (in moles) is:")
    # The final format is an ordered set, let's print it as a tuple
    result_tuple = tuple(moles_in_level)
    print(result_tuple)
    
    print("\n--- Breakdown ---")
    for i in range(len(n)):
        print(f"Level E{i+1}:")
        print(f"  Energy E_{i+1} = {E_n[i]:.4f} eV")
        print(f"  Degeneracy g_{i+1} = {g_n[i]}")
        print(f"  Probability P_{i+1} = {P_n[i]:.4f}")
        print(f"  Moles in Level = {moles_in_level[i]:.4f} mol")


# Execute the calculation and print the results
calculate_particle_distribution()

# The final answer in the required format
# To generate the answer string, we'll re-run the core logic
# since the function above prints more details.
epsilon_0_eV = 6.9e-3
T = 4200.0
k_B = 8.617333262145e-5
total_moles = 3.0
kT = k_B * T
n = np.arange(1, 6)
epsilon_n = epsilon_0_eV * (1 + n**n)
E_n = np.cumsum(epsilon_n)
g_n = 2 * n + 1
boltzmann_factors = g_n * np.exp(-E_n / kT)
Z1 = np.sum(boltzmann_factors)
P_n = boltzmann_factors / Z1
moles_in_level = total_moles * P_n
final_answer = tuple(moles_in_level)

# The <<<answer>>> format expects a clean, simple string representation of the tuple
print(f"\n<<<{final_answer}>>>")