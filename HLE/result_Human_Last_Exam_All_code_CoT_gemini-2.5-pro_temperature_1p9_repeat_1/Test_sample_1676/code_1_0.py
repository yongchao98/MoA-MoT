# First, you may need to install the RDKit library:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def identify_compound_3():
    """
    Identifies and characterizes the final product (Compound 3) of the synthesis.
    The final product is 1-isopropyl-4-methylcyclohex-3-ene-1-thiol.
    """
    # The SMILES string for 1-isopropyl-4-methylcyclohex-3-ene-1-thiol
    smiles_compound_3 = "CC1=CCC(C(C)C)(S)CC1"

    # Create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles_compound_3)

    if mol:
        # Calculate properties of the molecule
        formula = rdMolDescriptors.CalcMolFormula(mol)
        exact_mass = Descriptors.ExactMolWt(mol)

        print("--- Characterization of Compound 3 ---")
        print(f"Product Name: 1-isopropyl-4-methylcyclohex-3-ene-1-thiol")
        print(f"SMILES: {smiles_compound_3}")
        print(f"Chemical Formula: {formula}")
        print(f"Exact Molecular Weight: {exact_mass:.4f}")
    else:
        print("Error: Could not parse the SMILES string for Compound 3.")

if __name__ == "__main__":
    identify_compound_3()
