# Import necessary libraries from rdkit
from rdkit import Chem
from rdkit.Chem import Draw

def display_chemical_structures():
    """
    This function defines and displays the proposed structures for products A, B, and C.
    """
    # SMILES (Simplified Molecular Input Line Entry System) strings for the proposed structures
    smiles_A = "CC(=O)N1CC(N(C=C(C)C(=O)OC)C2)C2CCC1" # A proposed structure for Product A
    name_A = "Product A (C14H20N2O3)"

    smiles_B = "CC(=O)N1CCCC1(C(=O)O)c1cc[nH]c1" # N-acetyl-2-(pyrrol-2-yl)proline for Product B
    name_B = "Product B (C12H14N2O3)"
    
    smiles_C = "CC(=O)N1C2CCCN2C(=O)C2(N(C1=O)C2)" # A plausible tricyclic structure for Product C
    name_C = "Product C (C11H16N2O3)"

    # Create molecule objects from SMILES
    mol_A = Chem.MolFromSmiles(smiles_A)
    mol_B = Chem.MolFromSmiles(smiles_B)
    mol_C = Chem.MolFromSmiles(smiles_C)

    # Prepare for drawing
    mols = [mol_A, mol_B, mol_C]
    legends = [name_A, name_B, name_C]
    
    # Draw molecules to a PNG file and display it
    img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(300, 300), legends=legends)
    
    # In a real script, you would save this image or display it.
    # For this environment, we will print the SMILES strings as the definitive answer.
    print("Proposed Structures:")
    print("-" * 20)
    print(f"{name_A}:")
    print(f"SMILES: {smiles_A}\n")
    
    print(f"{name_B}:")
    print(f"SMILES: {smiles_B}\n")

    print(f"{name_C}:")
    print(f"SMILES: {smiles_C}\n")

# Execute the function
display_chemical_structures()

# To view the images if you run this code locally, you could use:
# img.save('products_ABC.png')
# from IPython.display import Image, display
# display(Image(filename='products_ABC.png'))