# First, we need to install the rdkit library if it's not already installed.
# You can do this by running: pip install rdkit-pypi
import sys
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("RDKit library not found. Please install it using: pip install rdkit-pypi")
    sys.exit(1)

# Step 1: Define the SMILES strings for the two products of the hydrolysis reaction.
# Based on the analysis, the products are Acetophenone and Diethylene Glycol.

acetophenone_smiles = "CC(=O)c1ccccc1"
diethyleneglycol_smiles = "OCCOCCO"

# Step 2: Create RDKit Mol objects from the SMILES strings.
mol_acetophenone = Chem.MolFromSmiles(acetophenone_smiles)
mol_diethyleneglycol = Chem.MolFromSmiles(diethyleneglycol_smiles)

# Step 3: Calculate the exact molecular weight (molar mass) for each molecule.
mass_acetophenone = Descriptors.MolWt(mol_acetophenone)
mass_diethyleneglycol = Descriptors.MolWt(mol_diethyleneglycol)

print(f"Product 1: Acetophenone")
print(f"SMILES: {acetophenone_smiles}")
print(f"Molar Mass: {mass_acetophenone:.4f} g/mol")
print("-" * 30)
print(f"Product 2: Diethylene Glycol")
print(f"SMILES: {diethyleneglycol_smiles}")
print(f"Molar Mass: {mass_diethyleneglycol:.4f} g/mol")
print("-" * 30)

# Step 4: Compare the molar masses and determine which is higher.
if mass_acetophenone > mass_diethyleneglycol:
    heavier_product_name = "Acetophenone"
    heavier_product_smiles = acetophenone_smiles
else:
    heavier_product_name = "Diethylene Glycol"
    heavier_product_smiles = diethyleneglycol_smiles

# Step 5: Print the final answer.
print(f"The product with the higher molar mass is {heavier_product_name}.")
print(f"Its SMILES string is:")
print(heavier_product_smiles)
