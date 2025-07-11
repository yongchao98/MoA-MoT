try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("RDKit library is not installed. Please install it using 'pip install rdkit-pypi'")
    # Provide dummy values to avoid crashing if rdkit is not present
    Descriptors = None

def solve_reaction():
    """
    Identifies the product with the higher molar mass from the hydrolysis of a ketal.

    The identity of the ketal is deduced from the fragments in the user-provided (but invalid) SMILES string.
    The assumed reaction is the acid hydrolysis of the ketal formed from acetophenone and diethylene glycol.
    The products are acetophenone and diethylene glycol.
    """
    if Descriptors is None:
        return

    # SMILES strings for the two hydrolysis products
    # Product 1: Acetophenone
    acetophenone_smiles = "CC(=O)c1ccccc1"
    # Product 2: Diethylene Glycol
    diethylene_glycol_smiles = "OCCOCCO"

    # Create RDKit molecule objects
    mol_acetophenone = Chem.MolFromSmiles(acetophenone_smiles)
    mol_diethylene_glycol = Chem.MolFromSmiles(diethylene_glycol_smiles)

    # Calculate molar masses
    mw_acetophenone = Descriptors.MolWt(mol_acetophenone)
    mw_diethylene_glycol = Descriptors.MolWt(mol_diethylene_glycol)
    
    # "Final equation" part: Printing the comparison
    print(f"Product 1: Acetophenone (SMILES: {acetophenone_smiles})")
    print(f"Molar Mass of Product 1: {mw_acetophenone:.4f}")
    
    print(f"\nProduct 2: Diethylene Glycol (SMILES: {diethylene_glycol_smiles})")
    print(f"Molar Mass of Product 2: {mw_diethylene_glycol:.4f}")

    # Determine which product has the higher molar mass
    if mw_acetophenone > mw_diethylene_glycol:
        heavier_product_smiles = acetophenone_smiles
    else:
        heavier_product_smiles = diethylene_glycol_smiles
    
    print(f"\nThe product with the higher molar mass is: {heavier_product_smiles}")

# Execute the function
solve_reaction()
