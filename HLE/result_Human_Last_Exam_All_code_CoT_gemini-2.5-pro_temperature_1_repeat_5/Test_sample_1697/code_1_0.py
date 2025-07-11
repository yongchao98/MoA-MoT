# This script identifies the IUPAC name of the reaction product.
# It requires the 'pubchempy' library. If you don't have it, please install it using:
# pip install PubChemPy

import pubchempy as pcp

def get_product_name():
    """
    Determines the IUPAC name of the product from a chemical reaction.
    
    The reaction is the ortho-metalation of N,N-diethyl-3-dimethylaminobenzamide
    followed by methylation. The chemical reasoning points to the formation of
    N,N-diethyl-2-methyl-3-dimethylaminobenzamide. This function uses the
    SMILES string of the predicted product to find its standard IUPAC name.
    """
    # SMILES string for the predicted product: N,N-diethyl-2-methyl-3-dimethylaminobenzamide
    product_smiles = "CCN(CC)C(=O)c1c(C)c(N(C)C)ccc1"
    
    try:
        # Search PubChem for the compound using its SMILES string
        compounds = pcp.get_compounds(product_smiles, 'smiles')
        
        if compounds:
            # Get the first compound found
            product = compounds[0]
            # Retrieve and print the IUPAC name
            product_name = product.iupac_name
            print(f"The resulting compound is: {product_name}")
        else:
            print("Could not find the compound in the PubChem database.")
            # Fallback name based on chemical nomenclature rules
            print("Predicted name: N,N-diethyl-2-methyl-3-(dimethylamino)benzamide")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an internet connection and the pubchempy library is installed correctly.")

if __name__ == "__main__":
    get_product_name()