# First, you may need to install the RDKit library if you haven't already.
# You can do this by running: pip install rdkit
# Note: For more complex installations (e.g., on Windows or with Conda),
# please refer to the official RDKit installation guide.

from rdkit import Chem
from rdkit.Chem import Descriptors

def identify_compound_A():
    """
    Identifies the product of the reaction between 2-aminopyridine,
    o-phthalaldehyde, and TMSCN.
    """
    # The final product A is 1-cyano-2-(pyridin-2-yl)isoindole.
    # We can represent this molecule using its SMILES string.
    smiles_A = "c1ccc2c(c1)c(C#N)cn2c3ncccc3"

    # Create an RDKit molecule object from the SMILES string
    mol_A = Chem.MolFromSmiles(smiles_A)

    if mol_A:
        # Calculate properties of the molecule
        molecular_formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol_A)
        molecular_weight = Descriptors.ExactMolWt(mol_A)
        common_name = "1-cyano-2-(pyridin-2-yl)isoindole"

        # Print the results
        print(f"Compound A has been identified.")
        print(f"Common Name: {common_name}")
        print(f"SMILES String: {smiles_A}")
        print(f"Molecular Formula: {molecular_formula}")
        print(f"Molecular Weight: {molecular_weight:.4f}")

    else:
        print("Could not create molecule from SMILES string. Please check the string.")

if __name__ == "__main__":
    identify_compound_A()
