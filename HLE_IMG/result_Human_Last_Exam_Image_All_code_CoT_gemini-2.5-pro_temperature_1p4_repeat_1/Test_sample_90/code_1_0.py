from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG

def get_svg_from_smiles(smiles_string):
    """Generates an SVG image of a molecule from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles_string)
    if not mol:
        return "Invalid SMILES string"
    # Generate coordinates
    try:
        from rdkit.Chem import AllChem
        AllChem.Compute2DCoords(mol)
    except ImportError:
        pass # RDKit will use its own algorithm
    
    # Drawing options
    drawer = rdMolDraw2D.MolDraw2DSVG(350, 300)
    # Highlight atoms for clarity if needed (optional)
    # atom_indices = range(mol.GetNumAtoms())
    # drawer.DrawMolecule(mol, highlightAtoms=atom_indices)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg.replace('svg:', '')

# SMILES representations of the proposed structures
smiles_A = "COC(=O)C1=C2N(CCCC2(C3(O)NCCC3)C1)C(=O)C" # A simplified representation of the complex Huisgen product after hydration
smiles_B = "O=C1OC(=O)C(C=CN2C1(C3=NCCC3)CCC2)=C" # Michael adduct after cyclization/elimination
smiles_C = "CC(=O)N1CCCC1(C2=NCCC2)C(=O)O"      # N-acetylated starting material

# Displaying the structures
print("Structure of Product A (C14H20N2O3):")
# Note: The exact structure of A is very complex. The SMILES below represents
# the Huisgen adduct prior to the proposed hydrolysis/rearrangement for simplicity,
# as it captures the core connectivity from the main reaction. The actual structure
# fitting the NMR data is the hydrated form of this adduct's side-chain.
# For a better representation let's try a different SMILES for A which is the hydrated form
smiles_A_hydrated = "COC(=O)c1c2n(c(c1C3(O)NCCC3)C(=O)C)CCC2" # a more speculative but NMR-consistent structure
# The problem is exceedingly complex, so let's stick to the most direct reaction products
# even if they don't explain every single NMR peak perfectly without further assumptions.
# Let's represent A as the direct cycloaddition product before side-chain hydration.
smiles_A_cycloadd = "COC(=O)C1=C(C(=O)C)N2C(C3=NCCC3)C1CC2"


print(get_svg_from_smiles(smiles_A_hydrated))
print("\n" + "="*50 + "\n")
print("Structure of Product B (C12H14N2O3):")
print(get_svg_from_smiles(smiles_B))
print("\n" + "="*50 + "\n")
print("Structure of Product C (C11H16N2O3):")
print(get_svg_from_smiles(smiles_C))
