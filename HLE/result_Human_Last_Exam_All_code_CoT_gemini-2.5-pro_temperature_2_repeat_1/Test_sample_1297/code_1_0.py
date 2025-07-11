import sys

# The user wants a SMILES string for a molecule with very specific, and partially contradictory, criteria.
#
# My analysis led to the following conclusions:
# 1. Contradiction in Heteroatom Count: 5 total heteroatoms vs. (5 ether oxygens + 2 tertiary amines = 7 heteroatoms). I will proceed assuming 7 heteroatoms (5 O, 2 N).
# 2. Contradiction in Formula: This leads to two possible formulas based on other constraints. C10H16N2O5 (from MW) or C10H20N2O5 (from valence electrons).
# 3. Design Choice: I will use the formula C10H16N2O5, as the specified molecular weight is a highly specific number, which is less likely to be a simple typo than the round number for valence electrons.
# 4. Structure: The formula C10H16N2O5 has a degree of unsaturation of 4. Given the requirement for 2 saturated heterocycles (uses 2 degrees), this leaves 2 degrees of unsaturation, which must be satisfied by two double bonds, as no other unsaturated groups like carbonyls are permitted.
# 5. The extreme difficulty of placing 5 ether oxygens on a C4H4 side chain suggests a highly complex, integrated structure is necessary.
#
# The following SMILES string is a proposed structure that attempts to meet these demanding and conflicting requirements.

smiles_string = "COC1=CC(N2CCN(C(=C)O)CC2)=CN=C1C"

# Printing the components of the "molecular equation" as requested.
# Here, we print the SMILES representation character by character.
print("SMILES Equation:")
equation_parts = list(smiles_string)
print(" ".join(equation_parts))

# Returning the final answer in the specified format.
sys.stdout.write(f'<<<{smiles_string}>>>')