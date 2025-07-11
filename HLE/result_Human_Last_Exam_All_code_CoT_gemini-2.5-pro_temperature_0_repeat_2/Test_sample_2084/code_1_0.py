# First, we need to install the RDKit library if it's not already installed.
# You can do this by running: pip install rdkit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("RDKit library not found. Please install it using 'pip install rdkit'")
    exit()

# The user provided an invalid SMILES string for the reactant.
# Based on the likely fragments (phenyl, methyl, ketal) and reaction conditions (acid hydrolysis),
# we assume the reaction is the hydrolysis of an acetophenone ketal.
# This reaction yields two products: acetophenone and a diol.

# SMILES string for Product 1: Acetophenone
smi_ketone = "O=C(C)c1ccccc1"

# SMILES string for Product 2: A simple diol. We'll use ethylene glycol as a representative example.
# The choice of a small diol is reasonable and won't change which product is heavier.
smi_diol = "OCCO"

# Create RDKit molecule objects from the SMILES strings
mol_ketone = Chem.MolFromSmiles(smi_ketone)
mol_diol = Chem.MolFromSmiles(smi_diol)

# Calculate the molar mass for each product
# This corresponds to the "numbers in the final equation" for comparison
mass_ketone = Descriptors.MolWt(mol_ketone)
mass_diol = Descriptors.MolWt(mol_diol)

print(f"Predicted Products of Hydrolysis:")
print(f"1. Acetophenone (SMILES: {smi_ketone})")
print(f"2. Ethylene Glycol (SMILES: {smi_diol})")
print("-" * 30)
print("Molar Mass Calculation:")
print(f"Molar Mass of Acetophenone: {mass_ketone:.2f} g/mol")
print(f"Molar Mass of Ethylene Glycol: {mass_diol:.2f} g/mol")
print("-" * 30)

# Compare the molar masses to find the heavier product
if mass_ketone > mass_diol:
    heavier_product_smi = smi_ketone
    heavier_product_name = "Acetophenone"
else:
    heavier_product_smi = smi_diol
    heavier_product_name = "Ethylene Glycol"

print(f"The product with the higher molar mass is: {heavier_product_name}")
print("SMILES string of the heavier product:")
print(heavier_product_smi)
