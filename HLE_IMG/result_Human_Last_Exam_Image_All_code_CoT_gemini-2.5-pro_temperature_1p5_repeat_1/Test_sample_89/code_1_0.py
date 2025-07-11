import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import io
import base64
from PIL import Image

def get_mol_image_b64(smiles):
    """Generates a base64 encoded PNG image of a molecule from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    AllChem.Compute2DCoords(mol)
    
    # Use IPython display to get a nice drawing
    drawer = Draw.rdMolDraw2D.MolDraw2DCairo(300, 200) # or MolDraw2DSVG
    drawer.drawOptions().addAtomIndices = False
    drawer.drawOptions().bondLineWidth = 2
    drawer.drawOptions().setAtomHighlightsAreCircles(True)
    drawer.drawOptions().setHighlightColour((0.8,0.8,0.8))

    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    
    # Get PNG data and encode it
    png_data = drawer.GetDrawingText()
    
    return base64.b64encode(png_data).decode('utf-8')


def print_final_answer():
    """
    This function determines and prints the structures for products A, B, and C
    based on the reaction descriptions provided by the user.
    """

    # --- Structure Determination based on analysis ---

    # Product A: Derived from the dihydropyrrole ring via an isocyanate intermediate,
    # followed by hydration, decarboxylation, and acetylation. This leads to an
    # open-chain N-acetylated amine.
    # SMILES string for N-(pent-4-en-1-yl)acetamide:
    smiles_A = "C=CCCCNC(C)=O"
    name_A = "N-(pent-4-en-1-yl)acetamide"
    formula_A = "C7H13NO"

    # Products B and C: The prompt states that both products are the "same tethered imide".
    # This imide is formed from the dihydropyrrole part reacting with acetic anhydride.
    # The most plausible structure is N-acetyl-2-pyrrolidinone, a cyclic N-acylated lactam.
    # SMILES string for N-acetyl-2-pyrrolidinone:
    smiles_B_C = "CC(=O)N1C(=O)CCC1"
    name_B_C = "N-acetyl-2-pyrrolidinone"
    formula_B_C = "C6H9NO2"


    # --- Printing the Final Answer ---
    print("The structures of the three products A, B, and C are as follows:\n")

    print("Product A:")
    print(f"Name: {name_A}")
    print(f"SMILES: {smiles_A}")
    print(f"Formula: {formula_A}\n")
    # In a real graphics-enabled environment, we would display the image.
    # For this text-based output, the name and SMILES are the structure representation.

    print("Product B:")
    print(f"Name: {name_B_C}")
    print(f"SMILES: {smiles_B_C}")
    print(f"Formula: {formula_B_C}\n")

    print("Product C:")
    print("Based on the reaction description, Product C is the same as Product B.")
    print(f"Name: {name_B_C}")
    print(f"SMILES: {smiles_B_C}")
    print(f"Formula: {formula_B_C}\n")

# Execute the function to print the solution.
print_final_answer()

# For a visual representation, here are the molecules encoded for display.
# (This part is illustrative and would require a capable frontend to render)
# image_A_b64 = get_mol_image_b64(smiles_A)
# image_B_C_b64 = get_mol_image_b64(smiles_B_C)
# print(f"Image data for A: data:image/png;base64,{image_A_b64}")
# print(f"Image data for B/C: data:image/png;base64,{image_B_C_b64}")
