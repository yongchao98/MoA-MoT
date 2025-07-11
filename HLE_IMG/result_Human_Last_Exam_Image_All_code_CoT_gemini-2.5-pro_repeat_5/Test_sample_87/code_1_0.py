# You may need to install RDKit: pip install rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors

def get_molecular_formula(smiles):
    """Calculates the molecular formula from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Invalid SMILES"
    mol = Chem.AddHs(mol)
    return Descriptors.rdkit_mol_formula(mol)

# SMILES strings for the three products
smiles_A = "COC(=O)C1=C(C(=O)C)C2C(N(C3=NCCCC3)CCC2)1"
smiles_B = "COC(=O)c1cn2c(c1)N(C3=NC(=O)CCC3)CCC2"
smiles_C = "O=C(O)C1N(C2=NC(C(=O)C)CCC2)CCC1"

# Target molecular formulas from the problem description
formula_A_target = "C14H20N2O3"
formula_B_target = "C12H14N2O3"
formula_C_target = "C11H16N2O3"

# Calculate formulas from proposed SMILES
formula_A_calculated = get_molecular_formula(smiles_A)
formula_B_calculated = get_molecular_formula(smiles_B)
formula_C_calculated = get_molecular_formula(smiles_C)

# Print the results
print("--- Product A ---")
print(f"Proposed Structure (SMILES): {smiles_A}")
print(f"Target Formula:   {formula_A_target}")
print(f"Calculated Formula: {formula_A_calculated}")
print(f"Match: {formula_A_target == formula_A_calculated}\n")

print("--- Product B ---")
print(f"Proposed Structure (SMILES): {smiles_B}")
print(f"Target Formula:   {formula_B_target}")
print(f"Calculated Formula: {formula_B_calculated}")
print(f"Match: {formula_B_target == formula_B_calculated}\n")

print("--- Product C ---")
print(f"Proposed Structure (SMILES): {smiles_C}")
print(f"Target Formula:   {formula_C_target}")
print(f"Calculated Formula: {formula_C_calculated}")
print(f"Match: {formula_C_target == formula_C_calculated}\n")