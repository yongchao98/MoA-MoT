import sys
import subprocess

# Ensure rdkit is installed
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("RDKit library not found. Installing it now...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "rdkit-pypi"])
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        print("RDKit installed successfully.")
    except Exception as e:
        print(f"Failed to install RDKit. Please install it manually using 'pip install rdkit-pypi'. Error: {e}")
        sys.exit(1)

# --- Plan Explanation ---
print("The chemical reaction is the acid-catalyzed hydrolysis of the provided ketal.")
print("The reaction breaks the starting material into two products: a ketone and a diol.")
print("Based on the fragments in the starting material's SMILES string, the two products are identified as Acetophenone and bis(hydroxymethyl)ether.")
print("\n--- Calculating Molar Masses ---")

# --- Product Definitions and Calculations ---
# Product 1: Acetophenone
smiles_product1 = "CC(=O)c1ccccc1"
name_product1 = "Acetophenone"
mol1 = Chem.MolFromSmiles(smiles_product1)
mass1 = Descriptors.ExactMolWt(mol1)

# Product 2: bis(hydroxymethyl)ether
smiles_product2 = "O(CO)2"
name_product2 = "bis(hydroxymethyl)ether"
mol2 = Chem.MolFromSmiles(smiles_product2)
mass2 = Descriptors.ExactMolWt(mol2)

# Print the molar masses used in the comparison
print(f"Product 1: {name_product1}")
print(f"   SMILES: {smiles_product1}")
print(f"   Molar Mass: {mass1:.4f} g/mol")

print(f"Product 2: {name_product2}")
print(f"   SMILES: {smiles_product2}")
print(f"   Molar Mass: {mass2:.4f} g/mol")

# --- Comparison and Final Answer ---
print("\n--- Comparison ---")
if mass1 > mass2:
    higher_mass_product_smiles = smiles_product1
    print(f"Acetophenone ({mass1:.4f}) has a higher molar mass than bis(hydroxymethyl)ether ({mass2:.4f}).")
else:
    higher_mass_product_smiles = smiles_product2
    print(f"bis(hydroxymethyl)ether ({mass2:.4f}) has a higher molar mass than Acetophenone ({mass1:.4f}).")

print("\nThe SMILES string for the product with the higher molar mass is:")
print(higher_mass_product_smiles)
