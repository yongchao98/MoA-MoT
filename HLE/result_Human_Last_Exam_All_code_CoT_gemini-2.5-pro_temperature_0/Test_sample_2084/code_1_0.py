# First, we need to install the RDKit library if it's not already installed.
# You can do this by running the following command in your terminal:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors

def solve_chemistry_problem():
    """
    This script solves the user's request by identifying the products of a chemical reaction,
    calculating their molar masses, and determining which is heavier.
    """
    # Step 1: Define the SMILES strings for the two likely products of the hydrolysis reaction.
    # Product 1: Acetophenone
    # Product 2: Diethylene Glycol
    product_info = {
        "Acetophenone": "CC(=O)c1ccccc1",
        "Diethylene Glycol": "OCCOCCO"
    }

    # Step 2: Calculate the molar mass for each product.
    product_masses = {}
    for name, smiles in product_info.items():
        molecule = Chem.MolFromSmiles(smiles)
        if molecule:
            molar_mass = Descriptors.MolWt(molecule)
            product_masses[name] = (molar_mass, smiles)

    # Step 3: Compare the molar masses to find the product with the higher mass.
    product1_name = "Acetophenone"
    product2_name = "Diethylene Glycol"
    
    mass1 = product_masses[product1_name][0]
    mass2 = product_masses[product2_name][0]

    if mass1 > mass2:
        heavier_product_name = product1_name
        heavier_product_mass = mass1
        heavier_product_smiles = product_masses[product1_name][1]
        lighter_product_mass = mass2
    else:
        heavier_product_name = product2_name
        heavier_product_mass = mass2
        heavier_product_smiles = product_masses[product2_name][1]
        lighter_product_mass = mass1

    # Step 4: Print the results, including the "equation" comparing the masses.
    print("The hydrolysis reaction produces Acetophenone and Diethylene Glycol.")
    print(f"Molar Mass of Acetophenone: {mass1:.2f} g/mol")
    print(f"Molar Mass of Diethylene Glycol: {mass2:.2f} g/mol")
    print("\nComparing the molar masses:")
    print(f"{heavier_product_mass:.2f} > {lighter_product_mass:.2f}")
    print(f"\nThe product with the higher molar mass is {heavier_product_name}.")
    print(f"The SMILES string for the heavier product is: {heavier_product_smiles}")

solve_chemistry_problem()