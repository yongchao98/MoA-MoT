from rdkit import Chem
from rdkit.Chem import Descriptors

def get_molecular_formula(smiles_string):
    """Calculates the molecular formula from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        return "Invalid SMILES"
    return Descriptors.rdMolDescriptors.CalcMolFormula(mol)

# Based on the analysis, we propose the following SMILES for the products.
# Note: The exact isomer (regio- and stereochemistry) can be ambiguous without more data,
# but the core connectivity and functional groups are based on the proposed reaction mechanism.

# Product A: C14H20N2O3 (Acetylated and reduced P)
# Structure: A pyrrolizidine core, substituted with R-group, acetyl, and methoxycarbonyl groups.
smiles_A = "COC(=O)C1C(C(=O)C)C2N(C3=NCCCC3)CCCC12"

# Product B: C12H14N2O3 (Oxidized P)
# Structure: A dihydropyrrolizine core, with R-group oxidized to a lactam.
smiles_B = "COC(=O)C1=CC2N(C3=NC(C(=O)CC3)=O)CCC=C12"

# Product C: C11H16N2O3 (Saponified and hydrated P)
# Structure: A dihydropyrrolizine core, with COOH group and hydrated R-group.
smiles_C = "O=C(O)C1=CC2N(C3NCCCC3O)CCC=C12"


print("Proposed Structures (SMILES notation) and their Molecular Formulas:\n")

# Product A
formula_A = get_molecular_formula(smiles_A)
print(f"Product A Structure (SMILES): {smiles_A}")
print(f"Expected Formula: C14H20N2O3")
print(f"Calculated Formula: {formula_A}\n")

# Product B
formula_B = get_molecular_formula(smiles_B)
print(f"Product B Structure (SMILES): {smiles_B}")
print(f"Expected Formula: C12H14N2O3")
print(f"Calculated Formula: {formula_B}\n")

# Product C
formula_C = get_molecular_formula(smiles_C)
print(f"Product C Structure (SMILES): {smiles_C}")
print(f"Expected Formula: C11H16N2O3")
print(f"Calculated Formula: {formula_C}\n")
