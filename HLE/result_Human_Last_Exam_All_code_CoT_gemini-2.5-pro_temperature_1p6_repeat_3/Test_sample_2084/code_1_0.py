# First, we need to install the RDKit library if it's not already installed.
# You can do this by running: pip install rdkit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("RDKit library not found.")
    print("Please install it using 'pip install rdkit-pypi'")
    exit()

# Step 1: Define the SMILES strings for the hydrolysis products.
# The malformed SMILES CC12COC(OC1)(OC2)C1=CC=CC=C1 suggests a spiroketal.
# Acid-catalyzed hydrolysis in water breaks it into a ketone and a diol.
# Based on the fragments, the most likely products are acetophenone and di(ethylene glycol).

# SMILES string for Acetophenone
ketone_smiles = 'CC(=O)c1ccccc1'
# SMILES string for Di(ethylene glycol)
diol_smiles = 'OCCOCCO'

# Step 2: Create molecule objects from the SMILES strings.
try:
    ketone_mol = Chem.MolFromSmiles(ketone_smiles)
    diol_mol = Chem.MolFromSmiles(diol_smiles)
    if ketone_mol is None or diol_mol is None:
        raise ValueError("Invalid SMILES string provided for products.")
except Exception as e:
    print(f"An error occurred while creating molecule objects: {e}")
    exit()

# Step 3: Calculate the molar mass (Molecular Weight) for each product.
ketone_mass = Descriptors.MolWt(ketone_mol)
diol_mass = Descriptors.MolWt(diol_mol)

# Step 4: Compare the molar masses and identify the product with the higher molar mass.
print("The hydrolysis reaction yields two products: Acetophenone and Di(ethylene glycol).")
print(f"1. Acetophenone (SMILES: {ketone_smiles})")
print(f"   Molar Mass = {ketone_mass:.4f} g/mol")
print(f"2. Di(ethylene glycol) (SMILES: {diol_smiles})")
print(f"   Molar Mass = {diol_mass:.4f} g/mol")
print("-" * 30)

if ketone_mass > diol_mass:
    higher_mass_product_name = "Acetophenone"
    higher_mass_product_smiles = ketone_smiles
else:
    higher_mass_product_name = "Di(ethylene glycol)"
    higher_mass_product_smiles = diol_smiles

print(f"The product with the higher molar mass is {higher_mass_product_name}.")
print("Its SMILES string is:")
print(higher_mass_product_smiles)
