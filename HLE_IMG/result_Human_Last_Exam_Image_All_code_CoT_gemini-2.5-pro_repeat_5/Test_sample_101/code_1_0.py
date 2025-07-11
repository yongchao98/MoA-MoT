# The rdkit library is required to run this code.
# You can install it via pip: pip install rdkit
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# --- Define the molecules involved in the reaction using SMILES strings ---

# Starting Material 1: 3-hydroxy-pyridine-2-carbaldehyde
reactant1_name = "3-hydroxy-pyridine-2-carbaldehyde"
reactant1_smiles = "O=Cc1c(O)cccn1"

# Starting Material 2: Aniline
reactant2_name = "Aniline"
reactant2_smiles = "c1ccccc1N"

# Final Product: Compound A is (3-hydroxypyridin-2-yl)(phenylamino)acetonitrile
# This is formed by the addition of HCN (from NaCN and a proton source) across the imine C=N bond.
product_A_name = "(3-hydroxypyridin-2-yl)(phenylamino)acetonitrile"
product_A_smiles = "N#CC(Nc1ccccc1)c2c(O)cccn2"

# --- Calculate Molecular Formulas ("the numbers in the equation") ---

# Create RDKit molecule objects from SMILES
mol_reactant1 = Chem.MolFromSmiles(reactant1_smiles)
mol_reactant2 = Chem.MolFromSmiles(reactant2_smiles)
mol_product_A = Chem.MolFromSmiles(product_A_smiles)

# Calculate the molecular formula for each molecule
# This function returns a string like "C13H11N3O"
formula_reactant1 = CalcMolFormula(mol_reactant1)
formula_reactant2 = CalcMolFormula(mol_reactant2)
formula_product_A = CalcMolFormula(mol_product_A)

# --- Print the final answer ---

print("The reaction produces Compound A, an alpha-aminonitrile.")
print(f"The identity of Compound A is: {product_A_name}\n")

print("Overall Reaction:")
print(f"{reactant1_name} + {reactant2_name} + NaCN  --->  Compound A\n")

print("--- Details of the final equation ---")
print("Molecule: 3-hydroxy-pyridine-2-carbaldehyde")
print(f"Molecular Formula: {formula_reactant1}\n")

print("Molecule: Aniline")
print(f"Molecular Formula: {formula_reactant2}\n")

print("Molecule: Compound A")
print(f"Molecular Formula: {formula_product_A}\n")

# The final answer is the structure of Compound A, which is (3-hydroxypyridin-2-yl)(phenylamino)acetonitrile.
# The structure can be represented by its SMILES string.
final_answer_smiles = product_A_smiles
# For the purpose of a simple answer format, we will provide the molecular formula.
final_answer_formula = formula_product_A
# print(f'<<<Answer: {final_answer_formula}>>>')