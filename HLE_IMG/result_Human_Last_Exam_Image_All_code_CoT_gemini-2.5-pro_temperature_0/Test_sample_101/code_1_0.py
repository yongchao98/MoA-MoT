# To run this code, you need to install the RDKit library:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def get_compound_properties(smiles_string, compound_name):
    """
    Calculates and prints properties of a chemical compound from its SMILES string.
    """
    try:
        # Create a molecule object from the SMILES string
        molecule = Chem.MolFromSmiles(smiles_string)
        if molecule is None:
            print(f"Error: Invalid SMILES string for {compound_name}.")
            return

        # Add hydrogens to the molecule model for accurate calculations
        molecule = Chem.AddHs(molecule)

        # Calculate properties
        molecular_formula = CalcMolFormula(molecule)
        exact_mass = Descriptors.ExactMolWt(molecule)

        # Print the results
        print(f"Analysis of Final Product: Compound A")
        print(f"Name: {compound_name}")
        print(f"SMILES: {smiles_string}")
        print(f"Molecular Formula (the final equation): {molecular_formula}")
        print(f"Exact Molecular Weight: {exact_mass:.4f}")

    except ImportError:
        print("RDKit library not found. Please install it using 'pip install rdkit'")
    except Exception as e:
        print(f"An error occurred: {e}")

# The SMILES string for the final product, (3-hydroxypyridin-2-yl)(phenylamino)acetonitrile
compound_A_smiles = "N#CC(Nc1ccccc1)c2ncccc2O"
compound_A_name = "(3-hydroxypyridin-2-yl)(phenylamino)acetonitrile"

get_compound_properties(compound_A_smiles, compound_A_name)
