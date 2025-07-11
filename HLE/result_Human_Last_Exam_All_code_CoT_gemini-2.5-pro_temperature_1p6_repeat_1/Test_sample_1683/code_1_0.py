from rdkit import Chem
from rdkit.Chem import Descriptors
import re

def identify_final_compound():
    """
    Identifies the final product (Compound 4) of the described chemical synthesis.
    The analysis concludes that Compound 4 is Benzophenone-2,2'-dicarboxylic acid.
    This script will define the molecule and print its key properties.
    """
    # The SMILES (Simplified Molecular Input Line Entry System) string for
    # Benzophenone-2,2'-dicarboxylic acid is O=C(O)c1ccccc1C(=O)c2ccccc2C(=O)O
    smiles_compound_4 = "O=C(O)c1ccccc1C(=O)c2ccccc2C(=O)O"
    
    # Create an RDKit molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles_compound_4)
    
    if mol:
        # IUPAC Name for this compound is 2-(2-carboxybenzoyl)benzoic acid
        iupac_name = "2-(2-carboxybenzoyl)benzoic acid"
        
        # Calculate the molecular formula
        mol_formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
        
        print("The final product, Compound 4, has the following properties:")
        print(f"IUPAC Name: {iupac_name}")
        print(f"Molecular Formula: {mol_formula}")
        
        # As requested, outputting the number for each atom in the final structure
        print("\nThe elemental composition of the final product is:")
        
        # Use regex to find all element-count pairs in the formula string
        # e.g., "C15H10O5" -> [('C', '15'), ('H', '10'), ('O', '5')]
        atomic_composition = re.findall(r'([A-Z][a-z]*)(\d*)', mol_formula)
        
        for element, count in atomic_composition:
            # If the count is not specified in the formula, it defaults to 1
            number = int(count) if count else 1
            print(f"Number of {element} atoms: {number}")

    else:
        print("Error: Could not generate the chemical structure.")

# Run the function to display the answer
identify_final_compound()
