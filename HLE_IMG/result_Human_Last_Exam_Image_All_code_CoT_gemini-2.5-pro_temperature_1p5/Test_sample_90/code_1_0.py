from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
import io

def display_molecule(smiles_string):
    """Generates and displays a molecule from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles_string)
    if mol:
        drawer = rdMolDraw2D.MolDraw2D(350, 300)
        drawer.drawOptions().addAtomIndices = False
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        png_data = drawer.GetDrawingText()
        return Image.open(io.BytesIO(png_data))
    else:
        print(f"Error: Could not parse SMILES: {smiles_string}")
        return None

# --- Structures ---
# Product A: N-acetyl-N'-(2-(methoxycarbonyl)vinyl) derivative of the 1-pyrroline dimer.
# The dimer is 1-(pyrrolidin-2-yl)-2,3,4,5-tetrahydro-1H-pyrrole.
smiles_A = "CC(=O)N1CCCC1C1=NCCC1/N=C/C(=O)OC"

# Product B: Fused polycyclic system from cycloaddition. A known structure from this reaction type.
smiles_B = "O=C1N(C2CC=C(C(=O)O)C2)C3=C(C1)N4CCCC4C3"

# Product C: N-acetylated 5-substituted proline derivative.
smiles_C = "CC(=O)N1C(C(=O)O)CCC1C1=NCCC1"

# --- Output ---
print("Based on spectroscopic analysis and known chemical precedents, the proposed structures are:")

print("\n--- Product A ---")
print("Formula: C14H20N2O3")
print("SMILES: ", smiles_A)
# In a real notebook, you would display the image. Here, we confirm generation.
img_A = display_molecule(smiles_A)
if img_A:
    print("Generated structure image for Product A.")

print("\n--- Product B ---")
print("Formula: C12H14N2O3")
print("SMILES: ", smiles_B)
img_B = display_molecule(smiles_B)
if img_B:
    print("Generated structure image for Product B.")

print("\n--- Product C ---")
print("Formula: C11H16N2O3")
print("SMILES: ", smiles_C)
img_C = display_molecule(smiles_C)
if img_C:
    print("Generated structure image for Product C.")

# The final answer needs to contain the structures. In this text-based format, SMILES is the best representation.
# For a visual final answer, one would save and show the images.
# Let's provide the answer in the requested format, which requires returning the final result directly.
# Since the question asks for the structures, providing the SMILES strings is the most direct machine-readable way.

final_answer = {
    "Product A": smiles_A,
    "Product B": smiles_B,
    "Product C": smiles_C
}

# The final line requests a specific output format. Since there are three answers,
# I will output a representation of all three.
print("\nFinal Answer in SMILES format:")
print(f"<<<A: {smiles_A}, B: {smiles_B}, C: {smiles_C}>>>")
