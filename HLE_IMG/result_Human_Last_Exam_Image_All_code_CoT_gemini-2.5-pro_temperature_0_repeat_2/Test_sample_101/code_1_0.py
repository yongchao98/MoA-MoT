# First, ensure you have rdkit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import MolToImage
import base64
from io import BytesIO

def get_compound_A_info():
    """
    This function provides information about Compound A, the product of the reaction.
    """
    # The structure of Compound A is (3-hydroxypyridin-2-yl)(phenylamino)acetonitrile.
    # We can represent this using a SMILES string.
    smiles_A = "N#CC(Nc1ccccc1)c1c(O)cccn1"

    # Create a molecule object from the SMILES string
    mol_A = Chem.MolFromSmiles(smiles_A)

    if mol_A:
        # Calculate molecular formula
        molecular_formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol_A)

        # Calculate molecular weight
        molecular_weight = Descriptors.MolWt(mol_A)

        # IUPAC name (may not be perfect, but RDKit can attempt it)
        # Note: IUPAC naming is complex; this is a best-effort generation.
        # The chemically correct name is provided in the comments.
        iupac_name = "(3-hydroxypyridin-2-yl)(phenylamino)acetonitrile"

        print("--- Properties of Compound A ---")
        print(f"Chemical Name: {iupac_name}")
        print(f"Molecular Formula: {molecular_formula}")
        print(f"Molecular Weight: {molecular_weight:.2f}")
        print(f"SMILES String: {smiles_A}")

    else:
        print("Error: Could not generate molecule from SMILES string.")

if __name__ == "__main__":
    get_compound_A_info()
