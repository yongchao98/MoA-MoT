# Import necessary libraries from RDKit for molecule handling and drawing.
# Please install it first if you haven't: pip install rdkit
from rdkit import Chem
from rdkit.Chem import Draw
import io
import base64
# This is a fallback for non-Jupyter environments to "display" an image.
from PIL import Image

def display_image_for_non_notebook(img):
    """
    A helper function to save the image to a file and notify the user,
    as direct display is not supported in all execution environments.
    """
    img.save("products_ABC.png")
    print("\n[Image containing the structures of Products A, B, and C has been saved as 'products_ABC.png']")


def identify_product_structures():
    """
    This function defines the structures of products A, B, and C using
    SMILES strings and generates a visual representation.
    """
    # SMILES strings for the proposed structures
    smiles_A = "CC(=O)N1CCCC1C1=CN(C=C(C)OC(=O))CCC1" # A simplified representation that fits the spirit of the derived structure
    smiles_A_plausible = "CC(=O)N1CCCC1C1CN(C=CC(=O)OC)C=C1" # N-acetyl-2-{1-[(E)-3-methoxy-3-oxoprop-1-en-1-yl]-2,5-dihydro-1H-pyrrol-2-yl}pyrrolidine
    smiles_B = "O=C1N2C(C(=O)OC)C=CC23CCCN3CC1" # methyl 2-oxo-2,3,5,6,7,8-hexahydro-1H-dipyrrolo[1,2-a:2',1'-c]pyrazine-1-carboxylate
    smiles_C = "CC(=O)N1CCCC1C(=O)O.C1=CCNC1" #This is a simplified way to represent the two parts for N-acetyl-2-(2,5-dihydropyrrol-2-yl)proline
    smiles_C_plausible = "CC(=O)N1[C@H](C(=O)O)CCC1[C@H]1C=CCN1" # N-acetyl-2-(2,5-dihydropyrrol-2-yl)pyrrolidine-2-carboxylic acid

    # Dictionary holding the product names and their SMILES strings
    structures = {
        "A": smiles_A_plausible,
        "B": smiles_B,
        "C": smiles_C_plausible
    }

    print("--- Proposed Structures of Products A, B, and C ---")

    mols = []
    legends = []
    for name, smi in structures.items():
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mols.append(mol)
            legends.append(f"Product {name}")
            print(f"\nProduct {name}:")
            print(f"Molecular Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
            print(f"SMILES: {smi}")
        else:
            print(f"\nCould not generate molecule from SMILES for Product {name}: {smi}")

    # Generate a grid image of the molecules
    if mols:
        img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(350, 350), legends=legends)
        display_image_for_non_notebook(img)

if __name__ == "__main__":
    identify_product_structures()
