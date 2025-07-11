import math

def solve_statistical_mechanics():
    """
    Solves the statistical mechanics problem to find the particle distribution.
    """
    # 1. Define Constants
    e0_meV = 6.9  # meV
    e0_eV = e0_meV * 1e-3  # Convert meV to eV
    T = 4200.0  # Kelvin
    k_B = 8.617333262e-5  # Boltzmann constant in eV/K
    total_moles = 3.0

    # Thermal energy
    k_B_T = k_B * T

    # 2. Calculate Energy Levels and Degeneracies
    num_levels = 5
    epsilons = []
    energies = []
    degeneracies = []

    cumulative_energy = 0
    for n in range(1, num_levels + 1):
        # Calculate epsilon_n
        epsilon_n = e0_eV * (1 + math.pow(n, n))
        epsilons.append(epsilon_n)
        
        # Calculate cumulative energy E_n
        cumulative_energy += epsilon_n
        energies.append(cumulative_energy)
        
        # Calculate degeneracy g_n
        g_n = 2 * n + 1
        degeneracies.append(g_n)

    # 3. Calculate the Single-Particle Partition Function (Z1)
    Z1 = 0
    Z1_terms = []
    for i in range(num_levels):
        g = degeneracies[i]
        E = energies[i]
        term = g * math.exp(-E / k_B_T)
        Z1_terms.append(term)
        Z1 += term

    print(f"Calculation Steps:")
    print(f"------------------")
    print(f"Constants:")
    print(f"  Base energy (e0) = {e0_eV:.6f} eV")
    print(f"  Temperature (T) = {T} K")
    print(f"  Boltzmann constant * Temperature (k_B*T) = {k_B_T:.6f} eV\n")
    
    print(f"Energy Levels (En) and Degeneracies (gn):")
    for i in range(num_levels):
        print(f"  Level n={i+1}: E{i+1} = {energies[i]:.6f} eV, g{i+1} = {degeneracies[i]}")
    print("\n")

    print(f"Single-Particle Partition Function (Z1):")
    print(f"  Z1 = sum(gn * exp(-En / (k_B*T))) for n=1 to 5")
    equation_str = "  Z1 = " + " + ".join([f"{t:.4f}" for t in Z1_terms])
    print(equation_str)
    print(f"  Z1 = {Z1:.6f}\n")

    # 4. Calculate Particle Distribution in Moles
    moles_distribution = []
    print("Most Probable Number of Moles in Each Energy Level (n_n):")
    print("  Formula: n_n = total_moles * (gn * exp(-En / (k_B*T))) / Z1\n")
    
    for i in range(num_levels):
        probability = Z1_terms[i] / Z1
        moles_n = total_moles * probability
        moles_distribution.append(moles_n)
        print(f"  Level E{i+1}:")
        print(f"    n{i+1} = {total_moles:.1f} * ({Z1_terms[i]:.4f} / {Z1:.4f})")
        print(f"    n{i+1} = {moles_n:.6f} moles\n")

    # 5. Output Final Result
    print("Final Answer:")
    print("The ordered set of the number of moles in each energy level (E1, E2, E3, E4, E5) is:")
    final_result = tuple(moles_distribution)
    print(f"{final_result}")

    # Hidden final answer tag
    print(f"\n<<<{final_result}>>>")

if __name__ == "__main__":
    solve_statistical_mechanics()