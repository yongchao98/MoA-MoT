# First, you might need to install the RDKit library:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def identify_compound_A():
    """
    This function identifies Compound A from the described reaction
    and prints its chemical properties.
    """
    # The structure of Compound A, 2-((cyano)(phenylamino)methyl)pyridin-3-ol,
    # can be represented by the following SMILES string.
    smiles_A = "c1c(C(C#N)Nc2ccccc2)c(O)ccn1"
    
    # Create an RDKit molecule object from the SMILES string
    mol_A = Chem.MolFromSmiles(smiles_A)
    
    if mol_A:
        # Calculate key molecular properties
        mol_formula = CalcMolFormula(mol_A)
        # We use ExactMolWt for a more precise mass
        exact_mw = Descriptors.ExactMolWt(mol_A)
        iupac_name = "2-((cyano)(phenylamino)methyl)pyridin-3-ol"

        print("--- Analysis of the Final Product (Compound A) ---")
        print(f"Reaction Pathway: This is a two-step synthesis involving an imine formation followed by a Strecker-type reaction.")
        print("\nCompound A has been identified as an alpha-aminonitrile with the following properties:")
        print(f"IUPAC Name: {iupac_name}")
        print(f"SMILES String: {smiles_A}")
        print(f"Molecular Formula: {mol_formula}")
        print(f"Exact Molecular Weight: {exact_mw:.4f}")
    else:
        print("Error: Could not generate the molecular structure from the provided SMILES string.")

if __name__ == "__main__":
    identify_compound_A()