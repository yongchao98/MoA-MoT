# First, ensure you have RDKit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def identify_compound_C():
    """
    This function identifies and characterizes the final product, Compound C.
    The structure is determined based on the provided reaction scheme.
    """
    
    # The structure of Compound C is 9-(diethylamino)-9-(2,4,6-trihydroxyphenyl)-9H-xanthene-1,3,6,8-tetraol.
    # We can represent this structure using its SMILES string.
    smiles_C = "CCN(CC)C1(c2c(O)cc(O)cc2Oc2cc(O)cc(O)c12)c1c(O)cc(O)cc1O"
    
    # Create an RDKit molecule object from the SMILES string
    mol_C = Chem.MolFromSmiles(smiles_C)
    
    if mol_C:
        # Calculate properties of Compound C
        molecular_formula = CalcMolFormula(mol_C)
        exact_molecular_weight = Descriptors.ExactMolWt(mol_C)
        
        # Print the results
        print("Compound C has been identified based on the reaction sequence.")
        print("-" * 50)
        print(f"Systematic Name: 9-(diethylamino)-9-(2,4,6-trihydroxyphenyl)-9H-xanthene-1,3,6,8-tetraol")
        print(f"SMILES String: {smiles_C}")
        print(f"Molecular Formula: {molecular_formula}")
        print(f"Exact Molecular Weight: {exact_molecular_weight:.4f}")
        print("-" * 50)
    else:
        print("Error: Could not generate the molecule from the SMILES string.")

if __name__ == '__main__':
    identify_compound_C()
