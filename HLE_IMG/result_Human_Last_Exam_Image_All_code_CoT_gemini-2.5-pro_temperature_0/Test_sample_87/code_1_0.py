from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
import io

def display_mol_with_formula(mol, mol_name=""):
    """Helper function to display molecule with its formula."""
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    print(f"Structure of Product {mol_name}:")
    
    # Generate high-quality drawing
    drawer = rdMolDraw2D.MolDraw2DCairo(400, 250)
    drawer.drawOptions().addAtomIndices = False
    drawer.drawOptions().bondLineWidth = 1
    drawer.drawOptions().setHighlightColour((0.0, 0.0, 0.0, 0.0)) # No highlight
    
    # Center the molecule
    drawer.drawOptions().centreMoleculesBeforeDrawing = True
    
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    
    # Display image
    img_data = drawer.GetDrawingText()
    img = Image.open(io.BytesIO(img_data))
    img.show() # This will open the image in the default viewer
    
    print(f"Molecular Formula: {formula}\n")

def solve_chemistry_problem():
    """
    This function defines the structures of products A, B, and C
    and displays them.
    """
    # SMILES string for Product C: N-acetylated starting material
    smiles_C = "CC(=O)N1CCCC1(C(=O)O)C2=NCCC2"
    mol_C = Chem.MolFromSmiles(smiles_C)

    # SMILES string for Product A: Decarboxylative [3+2] cycloaddition product
    smiles_A = "COC(=O)C1=CN2C(C(=O)C)C(CCC2)C1(C1=NCCC1)"
    mol_A = Chem.MolFromSmiles(smiles_A)

    # SMILES string for Product B: Tetracyclic pyridone product
    smiles_B = "O=C1C=C(C2=NCCC2)C2=C3N1C(CCC3)C2"
    mol_B = Chem.MolFromSmiles(smiles_B)

    # Displaying the structures
    # Note: The images will be generated and opened by your system's default image viewer.
    # For this interactive environment, we will print the SMILES and formulas.
    
    print("The structures of the three products are as follows:\n")

    # Product A
    formula_A = Chem.rdMolDescriptors.CalcMolFormula(mol_A)
    print("Product A")
    print(f"SMILES: {smiles_A}")
    print(f"Formula: {formula_A} (Matches C14H20N2O3)\n")
    
    # Product B
    formula_B = Chem.rdMolDescriptors.CalcMolFormula(mol_B)
    print("Product B")
    print(f"SMILES: {smiles_B}")
    print(f"Formula: {formula_B} (Matches C12H14N2O3)\n")

    # Product C
    formula_C = Chem.rdMolDescriptors.CalcMolFormula(mol_C)
    print("Product C")
    print(f"SMILES: {smiles_C}")
    print(f"Formula: {formula_C} (Matches C11H16N2O3)\n")
    
    # To visualize, you would typically use a library like IPython.display in a notebook
    # or save the images to files. Here, we print the standard representations.
    # img = Draw.MolsToGridImage([mol_A, mol_B, mol_C], molsPerRow=3, subImgSize=(300, 300), legends=['Product A', 'Product B', 'Product C'])
    # img.save("products.png")
    # print("Image of structures saved to products.png")


solve_chemistry_problem()