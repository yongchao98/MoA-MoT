# First, ensure you have rdkit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def get_molecular_formula(smiles: str) -> str:
    """Calculates the molecular formula from a SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "Invalid SMILES"
        return CalcMolFormula(mol)
    except Exception as e:
        return f"Error: {e}"

# --- Reactants ---
pivalaldehyde_name = "pivalaldehyde"
# SMILES for (CH3)3C-CHO
pivalaldehyde_smiles = "CC(C)(C)C=O"

ylide_name = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
# SMILES for (C6H5)3P=CH-CH2-(ortho-chlorophenyl)
ylide_smiles = "c1(Cl)ccccc1CC=[P](c2ccccc2)(c3ccccc3)c4ccccc4"

# --- Products ---
product_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
# SMILES for (CH3)3C-CH=CH-CH2-(ortho-chlorophenyl). Using E-isomer for representation.
product_smiles = "CC(C)(C)/C=C/Cc1c(Cl)cccc1"

byproduct_name = "triphenylphosphine oxide"
# SMILES for (C6H5)3P=O
byproduct_smiles = "O=[P](c1ccccc1)(c2ccccc2)c3ccccc3"

# Calculate molecular formulas
pivalaldehyde_formula = get_molecular_formula(pivalaldehyde_smiles)
ylide_formula = get_molecular_formula(ylide_smiles)
product_formula = get_molecular_formula(product_smiles)
byproduct_formula = get_molecular_formula(byproduct_smiles)

# Print the full reaction equation
print("The Wittig reaction is as follows:")
print("-" * 50)
print("Reactants:")
print(f"  1. {pivalaldehyde_name} (Formula: {pivalaldehyde_formula})")
print(f"  2. {ylide_name} (Formula: {ylide_formula})")
print("\nProducts:")
print(f"  1. {product_name} (Formula: {product_formula})")
print(f"  2. {byproduct_name} (Formula: {byproduct_formula})")
print("-" * 50)

print("\nFinal Reaction Equation:")
# The "numbers" in the equation are represented by the atomic counts in the molecular formulas.
print(f"{pivalaldehyde_name} ({pivalaldehyde_formula})  +  {ylide_name} ({ylide_formula})  -->  {product_name} ({product_formula})  +  {byproduct_name} ({byproduct_formula})")
