# The rdkit library is required for this script.
# If you don't have it, please install it via pip: pip install rdkit-pypi
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("Error: RDKit library not found. Please install it using 'pip install rdkit-pypi'")
    exit()

def solve_chemistry_problem():
    """
    Solves the problem by identifying the products of a hypothesized reaction,
    calculating their molar masses, and returning the SMILES string of the heavier product.
    """

    # Based on the analysis, the two products of the hydrolysis are Acetophenone and bis(4-hydroxybutyl)ether.
    product1_name = "Acetophenone"
    product1_smiles = "CC(=O)c1ccccc1"

    product2_name = "bis(4-hydroxybutyl)ether"
    product2_smiles = "OCCCCOCCCCO"

    # Create RDKit molecule objects
    mol1 = Chem.MolFromSmiles(product1_smiles)
    mol2 = Chem.MolFromSmiles(product2_smiles)

    # Calculate exact molar masses
    # This fulfills the request to "output each number in the final equation"
    mass1 = Descriptors.ExactMolWt(mol1)
    mass2 = Descriptors.ExactMolWt(mol2)

    print(f"Product 1: {product1_name}")
    print(f"SMILES: {product1_smiles}")
    print(f"Molar Mass: {mass1:.4f} g/mol")
    print("-" * 30)
    print(f"Product 2: {product2_name}")
    print(f"SMILES: {product2_smiles}")
    print(f"Molar Mass: {mass2:.4f} g/mol")
    print("-" * 30)
    
    # Determine the heavier product
    if mass1 > mass2:
        heavier_product_smiles = product1_smiles
        heavier_product_name = product1_name
    else:
        heavier_product_smiles = product2_smiles
        heavier_product_name = product2_name
    
    print(f"Conclusion: The product with the higher molar mass is {heavier_product_name}.")
    print("The SMILES string for this product is:")
    print(heavier_product_smiles)
    
    return heavier_product_smiles

# Execute the function and capture the final answer in the required format
final_smiles = solve_chemistry_problem()

print(f"\n<<<{final_smiles}>>>")
