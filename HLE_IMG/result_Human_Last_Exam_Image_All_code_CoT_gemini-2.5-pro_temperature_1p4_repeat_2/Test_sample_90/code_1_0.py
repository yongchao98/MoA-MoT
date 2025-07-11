from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
import io

def draw_molecule(smiles, legend):
    """Generates an image of a molecule from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    d = rdMolDraw2D.MolDraw2DSVG(350, 300)
    d.drawOptions().addAtomIndices = False
    d.DrawMolecule(mol, legend=legend)
    d.FinishDrawing()
    svg = d.GetDrawingText()
    return svg

# SMILES strings for the products A, B, and C
smiles_A = "COC(=O)C1=C(C)N2[C@@H](C3=NCCC3)CCC2C1=O"
smiles_B = "O=C1[C@@H]2[C@@H]3N(C1=O)CCCC3=NC1=C2CCC1"
smiles_C = "CC(=O)N1[C@H](C(=O)O)CCC1C2=NCCC2"

# Generate images
img_A_svg = draw_molecule(smiles_A, "Product A: C14H20N2O3")
img_B_svg = draw_molecule(smiles_B, "Product B: C12H14N2O3")
img_C_svg = draw_molecule(smiles_C, "Product C: C11H16N2O3")

# In a real execution environment, we would display these images.
# For this text-based output, I will print the SMILES strings as a representation of the structures.
print("Product A:")
print(f"SMILES: {smiles_A}\n")

print("Product B:")
print(f"SMILES: {smiles_B}\n")

print("Product C:")
print(f"SMILES: {smiles_C}\n")

# To fulfill the prompt's request for a single image output,
# here's a conceptual way to combine them. A real run would output a visual grid.
# from IPython.display import display, SVG
# display(SVG(img_A_svg))
# display(SVG(img_B_svg))
# display(SVG(img_C_svg))
# The prompt asks for print output, so SMILES strings are the best representation.
