# First, we need to install the rdkit library if it's not already installed.
# You can do this by running: pip install rdkit
# This script will proceed assuming rdkit is installed.

from rdkit import Chem
from rdkit.Chem import Descriptors

# Based on the analysis, the reaction is the acid hydrolysis of a ketal.
# We assume the starting molecule was a ketal formed from acetophenone and ethylene glycol.
# The two products of the hydrolysis are therefore acetophenone and ethylene glycol.

# SMILES string for Product 1: Acetophenone
product1_smiles = 'CC(=O)c1ccccc1'
# SMILES string for Product 2: Ethylene Glycol
product2_smiles = 'OCCO'

# Create RDKit molecule objects from the SMILES strings
mol_product1 = Chem.MolFromSmiles(product1_smiles)
mol_product2 = Chem.MolFromSmiles(product2_smiles)

# Calculate the molar mass for each product
# The numbers in the final output will be these molar masses.
mw_product1 = Descriptors.MolWt(mol_product1)
mw_product2 = Descriptors.MolWt(mol_product2)

# Determine which product has the higher molar mass
if mw_product1 > mw_product2:
    higher_mass_product_smiles = product1_smiles
    higher_mass_product_name = "Acetophenone"
else:
    higher_mass_product_smiles = product2_smiles
    higher_mass_product_name = "Ethylene Glycol"

# Print the details of the products and their molar masses for clarity
print(f"Product 1 is Acetophenone.")
print(f"SMILES: {product1_smiles}")
print(f"Molar Mass: {mw_product1:.4f} g/mol")
print("-" * 30)
print(f"Product 2 is Ethylene Glycol.")
print(f"SMILES: {product2_smiles}")
print(f"Molar Mass: {mw_product2:.4f} g/mol")
print("-" * 30)

# Print the final answer
print(f"The product with the higher molar mass is {higher_mass_product_name}.")
print("The SMILES string for this product is:")
print(higher_mass_product_smiles)