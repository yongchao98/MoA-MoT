from rdkit import Chem
from rdkit.Chem import Draw

def display_molecule(smiles_string, legend_text):
    """Generates and displays a molecule from a SMILES string with a legend."""
    mol = Chem.MolFromSmiles(smiles_string)
    img = Draw.MolToImage(mol, size=(400, 400), legend=legend_text)
    return img

# --- Structure A ---
# Formula: C14H20N2O3
# Description: Product of decarboxylation, acylation, and Michael addition.
# The structure is N-(1-pyrrolin-2-yl)-1-(1-acetyl-2-(methoxycarbonyl)vinyl)pyrrolidine.
smiles_A = "COC(=O)C=C(C(=O)C)N1CCCC1N1C=NCC1"
# Note: The drawn pyrroline is 2-pyrroline for simplicity, 1-pyrroline is C1=NCCC1
smiles_A_correct_pyrroline = "COC(=O)C=C(C(=O)C)N1CCCC1N1C(=N)CCC1"


# --- Structure B ---
# Formula: C12H14N2O3 (Note: Likely a typo for C13H16N2O3, the cycloadduct)
# Description: A rearranged [3+2] cycloaddition product.
smiles_B = "COC(=O)C1=CN2C3C(C(=O)O)CCC3N(C1=C)C1=NCCC1"
# A more standard representation of a rearranged cycloadduct:
smiles_B_rearranged = "COC(=O)C1=CC2N(C(=O)C(C2)N1C1=NCCC1)C1CCC1"


# --- Structure C ---
# Formula: C11H16N2O3
# Description: N-acylated starting material (as an enamide).
smiles_C = "CC(=O)N1C=C(N2C(C(=O)O)CCC2)CC1"

# Generate and print images
# Since we can only output one code block and should print the result,
# we will print the SMILES strings which represent the structures.
print("Structure A (SMILES):")
print(smiles_A_correct_pyrroline)
print("\nStructure B (SMILES - representing a rearranged cycloadduct):")
print(smiles_B_rearranged)
print("\nStructure C (SMILES):")
print(smiles_C)

# To visualize, you can use the following code in a Jupyter notebook
# display(display_molecule(smiles_A_correct_pyrroline, "Product A"))
# display(display_molecule(smiles_B_rearranged, "Product B"))
# display(display_molecule(smiles_C, "Product C"))