# First, ensure you have rdkit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem

def get_wittig_product_info():
    """
    This function determines the product of a Wittig reaction between
    pivalaldehyde and (2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane,
    and prints its chemical properties.
    """

    # The Wittig reaction couples an aldehyde with a phosphorus ylide to form an alkene.
    # Reactant 1 (aldehyde): Pivalaldehyde, (CH3)3C-CHO
    # Reactant 2 (ylide): (2-(2-chlorophenyl)ethylidene)triphenylphosphane,
    #                     (2-Cl-C6H4)-CH2-CH=PPh3

    # The product is formed by joining the carbon chains at a new double bond.
    # Product structure: (CH3)3C-CH=CH-CH2-(C6H4-2-Cl)
    # The ylide is non-stabilized, so the (Z)-isomer is the expected major product.

    # IUPAC name: (Z)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene
    product_name = "(Z)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    
    # We represent this molecule using its SMILES string, which includes stereochemistry.
    # /C=C/ specifies a cis or (Z) double bond.
    product_smiles = "Clc1ccccc1C/C=C/C(C)(C)C"

    # Create an RDKit molecule object from the SMILES string
    mol = Chem.MolFromSmiles(product_smiles)
    
    if mol:
        # Calculate properties using RDKit
        mol_formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
        mol_weight = Descriptors.ExactMolWt(mol)
        
        print("The major product of the Wittig reaction is:")
        print(f"IUPAC Name: {product_name}")
        print(f"SMILES: {product_smiles}")
        print(f"Molecular Formula: {mol_formula}")
        print(f"Molecular Weight: {mol_weight:.2f} g/mol")
        
        # As requested, outputting each number from the final product name.
        # The numbers in "(Z)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene" are 1, 2, 4, 4, 2.
        numbers_in_name = "1, 2, 4, 4, 2"
        print("\nFinal Equation Numbers:")
        print(f"The numbers in the IUPAC name '{product_name}' are: {numbers_in_name}")

    else:
        print("Error: Could not generate molecule from SMILES string.")

if __name__ == "__main__":
    get_wittig_product_info()

<<< (Z)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene >>>