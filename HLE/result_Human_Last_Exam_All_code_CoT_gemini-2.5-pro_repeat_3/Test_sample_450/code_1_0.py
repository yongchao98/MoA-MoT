import math

# Step 1: Define Constants
epsilon_0 = 6.9e-3  # eV
T = 4200.0          # Kelvin
k_B = 8.617333e-5   # Boltzmann constant in eV/K
N_moles_total = 3.0 # Total number of moles

# Calculate the thermal energy
kT = k_B * T

print(f"System Properties:")
print(f"epsilon_0 = {epsilon_0} eV")
print(f"T = {T} K")
print(f"Total Moles = {N_moles_total}")
print(f"Thermal Energy kT = {kT:.4f} eV\n")

# Step 2: Calculate Energy Levels and Degeneracies
energies = []
degeneracies = []
boltzmann_terms = []
current_E = 0.0

print("Calculating Energy Levels and Boltzmann Factors:")
for n in range(1, 6):
    # Calculate epsilon_n
    epsilon_n = epsilon_0 * (1 + n**n)
    
    # Calculate cumulative energy E_n
    current_E += epsilon_n
    energies.append(current_E)
    
    # Calculate degeneracy g_n
    g_n = 2 * n + 1
    degeneracies.append(g_n)
    
    # Calculate the term for the partition function sum
    term = g_n * math.exp(-current_E / kT)
    boltzmann_terms.append(term)
    
    print(f"--- Level n={n} ---")
    print(f"  Energy E_{n} = {current_E:.4f} eV")
    print(f"  Degeneracy g_{n} = {g_n}")
    print(f"  Boltzmann Term (g_n * exp(-E_n/kT)) = {term:.6f}")

# Step 3: Calculate the Single-Particle Partition Function Z_1
Z_1 = sum(boltzmann_terms)
print(f"\nSingle-particle partition function Z_1 = {Z_1:.4f}\n")

# Step 4: Calculate the Number of Moles in Each Level
moles_distribution = []
print("Calculating the most probable number of moles in each energy level:")
for i in range(5):
    n = i + 1
    # Probability P_n = term / Z_1
    # Moles_n = N_moles_total * P_n
    moles_n = N_moles_total * boltzmann_terms[i] / Z_1
    moles_distribution.append(moles_n)
    print(f"Moles in E_{n} = {N_moles_total:.1f} * ({boltzmann_terms[i]:.6f} / {Z_1:.4f}) = {moles_n:.4f} moles")

# Step 5: Output the Final Result
# The problem asks for the ordered set (E1, E2, E3, E4, E5) representing the number of moles.
# This is interpreted as the ordered set of moles for each energy level.
print("\nFinal ordered set representing the number of moles in each energy level (E1, E2, E3, E4, E5):")
# Create a tuple for the final answer format
final_answer_tuple = tuple(moles_distribution)
print(final_answer_tuple)

# Format for the final answer block
final_answer_str = f"({', '.join(f'{m:.4f}' for m in final_answer_tuple)})"
print(f"\n<<<{final_answer_str}>>>")