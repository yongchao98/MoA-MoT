# To run this script, you need to install the RDKit library.
# You can install it by running the following command in your terminal:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def identify_compound_4():
    """
    This function identifies and describes the final product (Compound 4) of the reaction sequence.
    Compound 4 is bis(2-methylphenyl)ketone.
    """
    # The chemical structure of bis(2-methylphenyl)ketone is represented by its SMILES string.
    smiles_string = "CC1=CC=CC=C1C(=O)C2=CC=CC=C2C"

    # Create a molecule object from the SMILES string
    molecule = Chem.MolFromSmiles(smiles_string)

    if molecule:
        # Calculate properties
        formula = rdMolDescriptors.CalcMolFormula(molecule)
        molecular_weight = Descriptors.MolWt(molecule)
        
        # IUPAC name can be complex to generate algorithmically for non-trivial molecules.
        # We will use the common and IUPAC-accepted name.
        iupac_name = "bis(2-methylphenyl)methanone"
        
        # Get atom counts for the final "equation"
        atom_counts = {}
        for atom in molecule.GetAtoms():
            symbol = atom.GetSymbol()
            atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
        
        # Print the results
        print("--- Analysis of Final Product (Compound 4) ---")
        print(f"The final product, Compound 4, is identified as: {iupac_name}")
        print(f"Also known as: di-o-tolyl ketone")
        print(f"\nSMILES String: {smiles_string}")
        print(f"Molecular Weight: {molecular_weight:.2f} g/mol")
        print(f"Molecular Formula: {formula}")
        
        print("\nFinal Equation (Atom Counts):")
        for atom, count in sorted(atom_counts.items()):
            print(f"Number of {atom} atoms: {count}")
    else:
        print("Error: Could not generate molecule from SMILES string.")

if __name__ == '__main__':
    identify_compound_4()
