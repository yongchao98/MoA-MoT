import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors

# The designed SMILES string based on the derived constraints.
smiles_string = "C1COCCN1COCOCN2COCOC2"

# The molecular formula derived from the SMILES string is C10H20N2O5.
# This formula meets the constraints of 17 heavy atoms, 100 valence electrons,
# 2 tertiary amines, and 5 ether oxygens.
# The following calculates the molecular weight for this formula, which is the only
# constraint that conflicts with the prompt's specifications.

# Monoisotopic atomic masses for the calculation.
mass_C = 12.00000
mass_H = 1.00783
mass_N = 14.00307
mass_O = 15.99491

# Atom counts in the designed molecule.
num_C = 10
num_H = 20
num_N = 2
num_O = 5

# Calculate the precise molecular weight.
molecular_weight = (num_C * mass_C) + (num_H * mass_H) + (num_N * mass_N) + (num_O * mass_O)

# The prompt asks to output each number in the final equation.
print("Molecular Weight Calculation:")
print(f"({num_C} * {mass_C}) + ({num_H} * {mass_H}) + ({num_N} * {mass_N}) + ({num_O} * {mass_O}) = {molecular_weight:.5f}")
print("\n---")

print("Designed Molecule SMILES:")
print(smiles_string)

# Verification using RDKit to confirm the structure meets the design goals.
# mol = Chem.MolFromSmiles(smiles_string)
# print(f"\nVerification:")
# print(f"Molecular Formula: {Descriptors.CalcMolFormula(mol)}")
# print(f"Heavy Atoms: {mol.GetNumHeavyAtoms()}")
# print(f"Valence Electrons: {Descriptors.NumValenceElectrons(mol)}")
# print(f"Number of Rings: {Descriptors.RingCount(mol)}")
# print(f"Rotatable Bonds: {Descriptors.NumRotatableBonds(mol)}")
# print(f"Hydrogen Bond Donors: {Descriptors.NumHDonors(mol)}")
# print(f"Hydrogen Bond Acceptors: {Descriptors.NumHAcceptors(mol)}")
# print(f"Ether Oxygens (SMARTS '[O]([#6])[#6]'): {len(mol.GetSubstructMatches(Chem.MolFromSmarts('[O]([#6])[#6]')))}")
# print(f"Tertiary Amines (SMARTS '[NX3]([#6])([#6])[#6]'): {len(mol.GetSubstructMatches(Chem.MolFromSmarts('[NX3]([#6])([#6])[#6]')))}")

<<<C1COCCN1COCOCN2COCOC2>>>