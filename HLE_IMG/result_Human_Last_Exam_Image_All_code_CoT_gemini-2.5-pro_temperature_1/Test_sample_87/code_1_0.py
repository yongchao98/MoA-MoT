# The following script uses the rdkit library.
# If you don't have it, please install it: pip install rdkit
# This script is for generating and verifying the chemical structures.
# The final output will be the structures in SMILES format.

from rdkit import Chem
from rdkit.Chem import Descriptors

def get_molecular_formula(smiles):
    """Calculates the molecular formula from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Invalid SMILES"
    return Descriptors.rdkit_mol_formula(mol)

# Proposed SMILES strings for the products
# Product A: Acylated [3+2] cycloadduct
smiles_A = "CC(=O)C1C(C(=O)OC)=CN(C2=NCCC2)[C@H]3CC[C@H]13"

# Product B: Oxidized and aromatized [3+2] cycloadduct
smiles_B = "COC(=O)c1cn(C2=[N+]([O-])CCC2)c2CCC12"

# Product C: N-Acetylated starting material (as a zwitterion)
smiles_C = "CC(=O)[N+]1=C(N2C(CCC2)C(=O)[O-])CCC1"

# Verify the formulas
formula_A_calculated = get_molecular_formula(smiles_A)
formula_B_calculated = get_molecular_formula(smiles_B)
formula_C_calculated = get_molecular_formula(smiles_C)

# Print the results
print("Based on the reaction mechanism and molecular formulas provided, the proposed structures are:")
print("-" * 30)

print("Product A:")
print(f"Proposed Molecular Formula: C14H20N2O3")
print(f"Calculated Formula from SMILES: {formula_A_calculated}")
print(f"SMILES representation: {smiles_A}")
print("Structure Description: The [3+2] cycloaddition product between the azomethine ylide (from decarboxylation of the starting material) and methyl propiolate, which is then C-acetylated on the resulting enamine ring.")
print("-" * 30)


print("Product B:")
print(f"Proposed Molecular Formula: C12H14N2O3")
print(f"Calculated Formula from SMILES: {formula_B_calculated}")
print(f"SMILES representation: {smiles_B}")
print("Structure Description: The primary [3+2] cycloaddition product that has undergone subsequent oxidation. This involves aromatization of the newly formed ring and oxidation of the substituent imine group to a nitrone.")
print("-" * 30)

print("Product C:")
print(f"Proposed Molecular Formula: C11H16N2O3")
print(f"Calculated Formula from SMILES: {formula_C_calculated}")
print(f"SMILES representation: {smiles_C}")
print("Structure Description: The starting material acetylated at the basic imine nitrogen.")
print("-" * 30)
