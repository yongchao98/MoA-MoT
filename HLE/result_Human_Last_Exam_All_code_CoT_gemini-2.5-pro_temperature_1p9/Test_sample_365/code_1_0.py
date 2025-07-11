import sys
import re

# This script requires the RDKit library to work.
# If you don't have it, you can install it using the command:
# pip install rdkit
try:
    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError:
    print("Error: The RDKit library is required but not found.")
    print("Please install it using: pip install rdkit")
    sys.exit(1)

def get_product_identity():
    """
    Analyzes the reaction and prints the identity of the final product.
    The reaction is an anionic oxy-Cope rearrangement, an isomerization process.
    """
    # The SMILES string for the starting material, (1S,2R,4S)-2-((S)-4-((tert-butyldimethylsilyl)oxy)cyclopent-1-en-1-yl)-7,7-dimethoxybicyclo[2.2.1]hept-5-en-2-ol.
    # A non-stereospecific SMILES is sufficient for calculating the molecular formula.
    start_material_smiles = "CC(C)(C)[Si](C)(C)OC1CC=C(C1)C2(O)C3C4C=CC3C4(OC)OC"

    # Create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(start_material_smiles)

    if mol is None:
        print("An error occurred while parsing the chemical structure.")
        return

    # Explain the nature of the product
    print("Product Analysis:")
    print("The reaction is an anionic oxy-Cope rearrangement, which is a type of isomerization.")
    print("The final product is a ketone with a rearranged bicyclo[8.3.0]tridecene skeleton.")
    print("-" * 30)

    # Calculate the molecular formula. Since it's an isomerization, the product has the same formula.
    product_formula = CalcMolFormula(mol)
    
    # Parse the formula string e.g., "C20H34O4Si" into a dictionary of element counts.
    # The regex finds pairs of (Element Symbol) and (Optional Count).
    atom_counts = {elem: int(num) if num else 1 for elem, num in re.findall(r'([A-Z][a-z]?)(\d*)', product_formula)}
    
    # Print the final product's chemical formula, detailing each element's count
    print("Product Chemical Formula:")
    print(f"Carbon (C): {atom_counts.get('C', 0)}")
    print(f"Hydrogen (H): {atom_counts.get('H', 0)}")
    print(f"Oxygen (O): {atom_counts.get('O', 0)}")
    print(f"Silicon (Si): {atom_counts.get('Si', 0)}")


if __name__ == "__main__":
    get_product_identity()
