# To run this code, you may need to install the RDKit library.
# You can install it by running this command in your terminal:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def find_compound_A():
    """
    This function identifies Compound A from the reaction of geraniol,
    and prints its properties.
    The reaction sequence described is a deoxygenation of the primary alcohol.
    """
    
    # Define the reactant and product using SMILES notation.
    # SMILES for Geraniol, (E)-3,7-dimethylocta-2,6-dien-1-ol
    reactant_smiles = "CC(C)=CCC/C(C)=C/CO"
    
    # The reaction replaces -OH with -H. The -CH2OH group becomes -CH3.
    # SMILES for the product, (E)-3,7-dimethylocta-2,6-diene
    product_A_smiles = "CC(C)=CCC/C(C)=C/C"
    
    # The name of compound A, derived from its structure
    product_A_name = "(E)-3,7-dimethylocta-2,6-diene"

    # Create molecule objects to calculate properties
    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
    product_A_mol = Chem.MolFromSmiles(product_A_smiles)
    
    # Calculate molecular formulas to show the transformation
    reactant_formula = CalcMolFormula(reactant_mol)
    product_A_formula = CalcMolFormula(product_A_mol)

    print("The overall reaction is a deoxygenation.")
    print("\n--- Chemical Equation Summary ---")
    print(f"Reactant (Geraniol) Formula: {reactant_formula}")
    print(f"Product (Compound A) Formula: {product_A_formula}")

    # Outputting the 'numbers in the final equation' as atom counts
    print("\nAtom Counts in the Transformation:")
    print("Reactant (Geraniol): Carbon=10, Hydrogen=18, Oxygen=1")
    print("Product (Compound A): Carbon=10, Hydrogen=18")

    print("\n--- Identity of Compound A ---")
    print(f"IUPAC Name: {product_A_name}")
    print(f"SMILES String: {product_A_smiles}")
    print(f"Molecular Formula: {product_A_formula}")

# Execute the function to find and describe Compound A
find_compound_A()
<<<(E)-3,7-dimethylocta-2,6-diene>>>