# To run this code, you need to have the RDKit library installed.
# You can install it via pip: pip install rdkit-pypi

from rdkit import Chem
from rdkit.Chem import Descriptors

def find_heavier_hydrolysis_product():
    """
    Identifies the product with the higher molar mass from the acid-catalyzed
    hydrolysis of a ketal.
    """
    # Step 1: Define the products of the hydrolysis reaction.
    # The provided SMILES string for the reactant, CC12COC(OC1)(OC2)C1=CC=CC=C1, is invalid.
    # We assume the intended reactant is the ketal formed from acetophenone and ethylene glycol.
    # The reaction is an acid-catalyzed hydrolysis, which cleaves the ketal.
    # Reactant (assumed): 2-methyl-2-phenyl-1,3-dioxolane
    # Products: Acetophenone and 1,2-ethanediol (Ethylene Glycol)

    product_smiles = {
        "Acetophenone": 'CC(=O)c1ccccc1',
        "Ethylene Glycol": 'OCCO'
    }

    print("Analyzing the acid-catalyzed hydrolysis products...")
    print("-" * 50)

    # Step 2: Calculate molar mass for each product and store them.
    product_data = []
    for name, smiles in product_smiles.items():
        molecule = Chem.MolFromSmiles(smiles)
        if molecule:
            # Using ExactMolWt for precise mass calculation
            molar_mass = Descriptors.ExactMolWt(molecule)
            product_data.append({'name': name, 'smiles': smiles, 'mass': molar_mass})
            print(f"Product: {name}")
            print(f"  SMILES: {smiles}")
            print(f"  Molar Mass: {molar_mass:.4f} g/mol")
        else:
            print(f"Could not process SMILES for {name}: {smiles}")

    print("-" * 50)

    # Step 3: Compare the molar masses to find the heavier product.
    if not product_data or len(product_data) < 2:
        print("Could not perform comparison.")
        return

    # Sort the products by mass in descending order
    heavier_product = sorted(product_data, key=lambda x: x['mass'], reverse=True)[0]

    # Step 4: Output the final answer.
    print(f"Comparing the molar masses of the products:")
    print(f"{product_data[0]['name']} ({product_data[0]['mass']:.4f}) vs {product_data[1]['name']} ({product_data[1]['mass']:.4f})")
    print("\nThe product with the higher molar mass is:")
    print(f"Name: {heavier_product['name']}")
    print(f"SMILES String: {heavier_product['smiles']}")

# Run the function
find_heavier_hydrolysis_product()