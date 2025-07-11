from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
import io

def display_molecule_with_info(smiles, label, formula, description):
    """Generates and displays an image of a molecule with labels."""
    mol = Chem.MolFromSmiles(smiles)
    
    # Prepare the molecule for drawing
    drawer = rdMolDraw2D.MolDraw2DCairo(500, 400)
    drawer.drawOptions().addAtomIndices = False
    drawer.drawOptions().bondLineWidth = 2
    drawer.drawOptions().setAtomBlack()
    
    # Set title
    title = f"Product {label}: {formula}\n{description}"
    drawer.drawOptions().legendFontSize=20

    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, legend=title)
    drawer.FinishDrawing()
    
    # Get image as bytes and display
    bio = io.BytesIO(drawer.GetDrawingText())
    img = Image.open(bio)
    
    # In a real script, you would save or show the image
    # For this environment, we will save to a file and print the filename.
    filename = f"product_{label}.png"
    img.save(filename)
    print(f"Structure of Product {label} ({formula}) has been saved to '{filename}'")
    # For inline display in notebooks, one might use: display(img)

def solve_chemistry_problem():
    """
    Solves the chemical structure problem by defining the SMILES strings for
    products A, B, and C and generating images.
    """
    # Product C: N-acetylated starting material
    # Formula: C11H16N2O3
    smiles_C = "CC(=O)N1CCC[C@]1(C(=O)O)C2=NCCC2"
    desc_C = "N-acetyl-2-(1'-pyrrolin-2'-yl)pyrrolidine-2-carboxylic acid"
    
    # Product A: Tricyclic product from [3+2] cycloaddition of the N-acetylated ylide and methyl propiolate
    # Formula: C14H20N2O3
    smiles_A = "COC(=O)C1=C[C@@H]2[C@H](N(C(=O)C)C1)CC[C@]2(C1=NCCC1)"
    desc_A = "Cycloadduct from N-acetylated ylide and methyl propiolate"

    # Product B: Pyrrolizinone derivative from cycloaddition, oxidation, and lactam formation
    # Formula: C12H14N2O3
    smiles_B = "COC(=O)C1=CC(=O)C2=C(C=N2)C1N1CCCC1=O"
    desc_B = "Pyrrolizinone substituted with a pyrrolidone lactam"

    print("--- Proposed Structures for the Reaction Products ---")
    
    # Generate images for each product
    display_molecule_with_info(smiles_C, "C", "C₁₁H₁₆N₂O₃", desc_C)
    display_molecule_with_info(smiles_A, "A", "C₁₄H₂₀N₂O₃", desc_A)
    display_molecule_with_info(smiles_B, "B", "C₁₂H₁₄N₂O₃", desc_B)
    
    print("\nNote: The image files have been generated in your current directory.")
    
# Execute the solver
solve_chemistry_problem()
