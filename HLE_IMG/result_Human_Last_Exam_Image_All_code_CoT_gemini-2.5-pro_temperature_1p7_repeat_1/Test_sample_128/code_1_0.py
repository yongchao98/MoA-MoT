# To run this code, you may need to install the RDKit library.
# You can do this by running: pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# Step 1: Define the molecules involved in the reaction using SMILES strings.
# The reaction is a 1:1:1 condensation of the three main components.
smiles_amine = "Nc1ncccc1"          # Reactant 1: 2-aminopyridine
smiles_aldehyde = "O=Cc1ccccc1C=O"  # Reactant 2: o-phthalaldehyde
smiles_cyanide_source = "[H]C#N"    # Reactant 3: Hydrogen Cyanide (from TMSCN and water)

# Step 2: Define the product A and byproduct based on the reaction mechanism.
# Product A is 3-hydroxy-2-(pyridin-2-yl)isoindoline-1-carbonitrile.
smiles_product_A = "c1ccc2c(c1)C(C#N)N(c3ncccc3)C(O)2"
smiles_byproduct = "O"              # Byproduct: Water

# Step 3: Create RDKit molecule objects from the SMILES strings.
mol_amine = Chem.MolFromSmiles(smiles_amine)
mol_aldehyde = Chem.MolFromSmiles(smiles_aldehyde)
mol_cyanide = Chem.MolFromSmiles(smiles_cyanide_source)
mol_product_A = Chem.MolFromSmiles(smiles_product_A)
mol_water = Chem.MolFromSmiles(smiles_byproduct)

# Step 4: Calculate the molecular formula for each component.
formula_amine = CalcMolFormula(mol_amine)
formula_aldehyde = CalcMolFormula(mol_aldehyde)
formula_cyanide = CalcMolFormula(mol_cyanide)
formula_product_A = CalcMolFormula(mol_product_A)
formula_water = CalcMolFormula(mol_water)

# Step 5: Display the balanced chemical equation with coefficients and molecular formulas.
print("The reaction is a three-component condensation with the elimination of one water molecule.")
print("\n--- Balanced Chemical Equation ---")
# The "numbers in the final equation" are the coefficients and the molecular formulas.
print(f"1 {formula_amine}  +  1 {formula_aldehyde}  +  1 {formula_cyanide}   ->   1 {formula_product_A}  +  1 {formula_water}")
print("-" * 35)

# Step 6: Identify and describe Compound A.
product_name = "3-hydroxy-2-(pyridin-2-yl)isoindoline-1-carbonitrile"
print(f"\nCompound A is identified as: {product_name}")
print(f"Molecular Formula: {formula_product_A}")

exact_mass_A = Descriptors.ExactMolWt(mol_product_A)
print(f"Exact Molecular Weight: {exact_mass_A:.4f}")
