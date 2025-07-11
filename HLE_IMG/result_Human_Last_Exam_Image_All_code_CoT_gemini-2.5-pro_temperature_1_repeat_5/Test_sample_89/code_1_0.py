import rdkit
from rdkit import Chem
from rdkit.Chem import Draw

def get_product_structures():
    """
    This function defines the SMILES strings for products A, B, and C based on the reaction description
    and prints them to the console.
    """
    # Product A: N-(3,4-dihydro-2H-pyrrol-5-yl)acetamide
    # Formed from the isocyanate pathway, followed by hydration, decarboxylation, and acetylation.
    smiles_A = "CC(=O)NC1=NCCC1"

    # Product B: N-(3,4-dihydro-2H-pyrrol-5-yl)maleimide
    # Formed from the Michael addition of propiolate, which then forms a maleimide ring.
    smiles_B = "O=C1C=CC(=O)N1C2=NCCC2"

    # Product C: N-(3,4-dihydro-2H-pyrrol-5-yl)diacetamide
    # Formed from the reaction with acetic anhydride, leading to diacetylation of the amine fragment.
    smiles_C = "CC(=O)N(C(=O)C)C1=NCCC1"

    print("The chemical structures of the products are represented by their SMILES strings:")
    print("-" * 60)
    print(f"Product A: {smiles_A}")
    print(f"Product B: {smiles_B}")
    print(f"Product C: {smiles_C}")
    print("-" * 60)

    # Note: To visualize these structures, you can use online tools that convert SMILES to images,
    # or use a cheminformatics library like RDKit as shown below (code is commented out).
    #
    # mol_A = Chem.MolFromSmiles(smiles_A)
    # mol_B = Chem.MolFromSmiles(smiles_B)
    # mol_C = Chem.MolFromSmiles(smiles_C)
    # img = Draw.MolsToGridImage([mol_A, mol_B, mol_C], molsPerRow=3, legends=['Product A', 'Product B', 'Product C'])
    # img.save('products.png') # This would save the image to a file.

if __name__ == '__main__':
    get_product_structures()
