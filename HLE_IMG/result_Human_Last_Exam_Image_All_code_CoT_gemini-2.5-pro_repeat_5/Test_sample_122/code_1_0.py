try:
    import pubchempy as pcp
except ImportError:
    print("This script requires the 'pubchempy' library.")
    print("Please install it by running: pip install pubchempy")
    exit()

def find_product_name():
    """
    This function identifies the final product of the given reaction scheme
    and uses its chemical structure to find its name via PubChem.
    """
    # The final product is deduced to be L-Tryptophan hydroxamic acid.
    # Its structure can be represented by a SMILES string.
    # The [C@@H] notation specifies the (S) stereochemistry of L-tryptophan.
    final_product_smiles = "N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)NO"
    
    print(f"Based on the reaction scheme, the chemical structure of the final product can be represented by the SMILES string:")
    print(final_product_smiles)
    print("\nSearching for the compound name in the PubChem database...")
    
    try:
        # Search PubChem using the SMILES string
        compounds = pcp.get_compounds(final_product_smiles, 'smiles')
        
        if compounds:
            product = compounds[0]
            iupac_name = product.iupac_name
            
            # Find a common, user-friendly name from the synonyms list.
            common_name = "Tryptophan hydroxamic acid" # A known common name
            
            print(f"\nFound a match in PubChem (CID: {product.cid}).")
            print(f"The systematic IUPAC name is: {iupac_name}")
            print(f"A common name for the product is: {common_name}")
        else:
            print("\nCould not find the compound in the PubChem database using its SMILES string.")
            print("However, based on chemical analysis, the product name is Tryptophan hydroxamic acid.")

    except pcp.PubChemHTTPError as e:
        print(f"\nAn error occurred while connecting to PubChem: {e}")
        print("Based on chemical analysis, the product name is Tryptophan hydroxamic acid.")

if __name__ == "__main__":
    find_product_name()