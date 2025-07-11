# To run this code, you may need to install the pubchempy library first:
# pip install pubchempy

import pubchempy as pcp

def identify_product():
    """
    Identifies the product from the reaction scheme and prints its details from PubChem.
    The reaction sequence converts the carboxylic acid of Tryptophan into a hydroxamic acid.
    The final product is Tryptophan hydroxamic acid.
    """
    product_name = "Tryptophan hydroxamic acid"
    
    try:
        # Search for the compound by its common name
        compounds = pcp.get_compounds(product_name, 'name')
        if not compounds:
            print(f"Could not find information for the product: {product_name}")
            return

        product = compounds[0]
        
        print(f"The final product of the synthesis is: {product_name}")
        print("\n--- Chemical Details ---")
        print(f"PubChem CID: {product.cid}")
        
        # The stereochemistry of the starting material is preserved. The IUPAC name reflects this.
        # The IUPAC name is (2S)-2-amino-N-hydroxy-3-(1H-indol-3-yl)propanamide
        print(f"IUPAC Name: {product.iupac_name}")
        
        print(f"Molecular Formula: {product.molecular_formula}")
        print(f"Molecular Weight: {product.molecular_weight}")

    except Exception as e:
        print(f"An error occurred while fetching data from PubChem: {e}")
        print("\nBased on manual analysis of the reaction scheme:")
        print(f"The product name is: {product_name}")
        print("IUPAC Name: (2S)-2-amino-N-hydroxy-3-(1H-indol-3-yl)propanamide")

if __name__ == "__main__":
    identify_product()
