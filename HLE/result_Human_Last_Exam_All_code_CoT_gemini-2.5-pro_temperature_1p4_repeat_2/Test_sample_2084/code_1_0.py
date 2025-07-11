# The rdkit library is a powerful toolkit for cheminformatics.
# If you don't have it installed, you can do so via pip:
# pip install rdkit-pypi

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    # This block provides a fallback if rdkit is not installed.
    print("Warning: RDKit library not found. Calculations will proceed with pre-computed values.")
    print("To run this dynamically, please install RDKit: pip install rdkit-pypi\n")
    
    product1_name = "Acetophenone"
    product1_smiles = "CC(=O)c1ccccc1"
    molar_mass_product1 = 120.15

    product2_name = "Ethylene Glycol"
    product2_smiles = "OCCO"
    molar_mass_product2 = 62.07
    
    print("Based on chemical analysis, the reaction is a ketal hydrolysis.")
    print(f"The likely products are {product1_name} and a simple diol like {product2_name}.")
    print("\nComparing their molar masses:")
    print(f"Molar mass of {product1_name} ({product1_smiles}): {molar_mass_product1:.2f} g/mol")
    print(f"Molar mass of {product2_name} ({product2_smiles}): {molar_mass_product2:.2f} g/mol")

    heavier_product_smiles = product1_smiles

    print(f"\nSince {molar_mass_product1:.2f} g/mol > {molar_mass_product2:.2f} g/mol, the heavier product is {product1_name}.")
    print("\nThe SMILES string for this product is:")
    print(heavier_product_smiles)

else:
    # This block executes if RDKit is found and will perform live calculations.
    # Step 1: Define the potential products based on chemical analysis.
    # The parent ketone is identified as acetophenone, and the simplest diol is ethylene glycol.
    product1_name = "Acetophenone"
    product1_smiles = "CC(=O)c1ccccc1"

    product2_name = "Ethylene Glycol"
    product2_smiles = "OCCO"

    # Step 2: Create RDKit molecule objects from their SMILES strings.
    mol_product1 = Chem.MolFromSmiles(product1_smiles)
    mol_product2 = Chem.MolFromSmiles(product2_smiles)

    # Step 3: Calculate the exact molar weight for each molecule.
    molar_mass_product1 = Descriptors.MolWt(mol_product1)
    molar_mass_product2 = Descriptors.MolWt(mol_product2)

    # Step 4: Compare molar masses and print the results as requested.
    print("The reaction is an acid-catalyzed hydrolysis of a ketal.")
    print(f"The most plausible products are {product1_name} and {product2_name}.")
    print("\nComparing the molar masses of the two products:")
    
    # Outputting the numbers for the comparison as requested.
    print(f"Molar Mass of {product1_name} ({product1_smiles}): {molar_mass_product1:.4f} g/mol")
    print(f"Molar Mass of {product2_name} ({product2_smiles}): {molar_mass_product2:.4f} g/mol")

    if molar_mass_product1 > molar_mass_product2:
        heavier_product_name = product1_name
        heavier_product_smiles = product1_smiles
    else:
        heavier_product_name = product2_name
        heavier_product_smiles = product2_smiles

    print(f"\nThe product with the higher molar mass is {heavier_product_name}.")
    print("Its SMILES string is:")
    print(heavier_product_smiles)
