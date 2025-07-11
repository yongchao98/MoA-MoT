# First, ensure you have the pubchempy library installed.
# You can install it using pip:
# pip install pubchempy

import pubchempy as pcp

def find_product_name():
    """
    Identifies the IUPAC name of the final product from the reaction sequence.
    The final product is L-Tryptophan hydroxamate, as deduced from the reaction steps.
    """
    # The canonical SMILES string for L-Tryptophan hydroxamate
    # (S)-2-amino-N-hydroxy-3-(1H-indol-3-yl)propanamide
    product_smiles = "C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)NO)N"

    try:
        # Search PubChem for the compound using its SMILES string
        compounds = pcp.get_compounds(product_smiles, 'smiles')

        if compounds:
            # Get the first matching compound
            product = compounds[0]
            # Retrieve and print the IUPAC name
            iupac_name = product.iupac_name
            if iupac_name:
                print(f"The name of the product is: {iupac_name}")
            else:
                print("IUPAC name not found, a common synonym is: L-Tryptophan hydroxamate")
        else:
            # Fallback if the compound is not found via SMILES
            print("Could not find the product on PubChem. Based on chemical analysis, the name is (2S)-2-amino-N-hydroxy-3-(1H-indol-3-yl)propanamide.")

    except Exception as e:
        print(f"An error occurred while connecting to PubChem: {e}")
        print("Based on chemical analysis, the name of the product is (2S)-2-amino-N-hydroxy-3-(1H-indol-3-yl)propanamide.")

if __name__ == "__main__":
    find_product_name()