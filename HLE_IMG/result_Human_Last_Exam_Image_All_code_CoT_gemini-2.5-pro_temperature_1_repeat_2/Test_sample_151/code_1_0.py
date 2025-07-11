# To run this code, you may need to install the pubchempy library first:
# pip install pubchempy

import pubchempy as pcp

def identify_product_name():
    """
    Identifies the IUPAC name of the final product from the reaction sequence.
    The final product is determined to be ethyl 2,5-dihydrothiophene-3-carboxylate,
    which has the SMILES string "CCOC(=O)C1=CSCC1".
    This script uses pubchempy to retrieve the IUPAC name from this SMILES string.
    """
    product_smiles = "CCOC(=O)C1=CSCC1"

    print(f"Analyzing the product with SMILES: {product_smiles}\n")

    try:
        # Fetch compound information from PubChem using its SMILES string
        compounds = pcp.get_compounds(product_smiles, 'smiles')

        if compounds:
            product = compounds[0]
            iupac_name = product.iupac_name

            print(f"The IUPAC name of the final product is: {iupac_name}")

            # Explain the numbers in the name as requested
            print("\n--- Name Breakdown ---")
            print("Parent Heterocycle: The core is a five-membered ring with a sulfur atom ('thiophene').")
            print("Saturation (positions 2, 5): The ring is partially hydrogenated, with saturated carbons at positions 2 and 5, hence '2,5-dihydro'.")
            print("Substituent (position 3): The main functional group, the ester, is at position 3, hence '-3-carboxylate'.")
            print("Ester Group: The ester is an ethyl ester, hence 'ethyl'.")

        else:
            print("Could not retrieve compound information from PubChem.")
            print("Based on chemical analysis, the name is: ethyl 2,5-dihydrothiophene-3-carboxylate")

    except Exception as e:
        print(f"An error occurred while contacting PubChem: {e}")
        print("Based on chemical analysis, the name is: ethyl 2,5-dihydrothiophene-3-carboxylate")

if __name__ == "__main__":
    identify_product_name()
