# First, ensure you have the pubchempy library installed:
# pip install pubchempy

import pubchempy as pcp

def find_chemical_name():
    """
    This function determines the name of the final product of the given reaction sequence.
    The final product is Tryptophan hydroxamate.
    Its structure can be represented by the SMILES string: N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)NO
    We will use this SMILES string to query the PubChem database for its common name.
    """
    try:
        # SMILES string for L-Tryptophan hydroxamate
        smiles = 'N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)NO'
        
        # Search PubChem for the compound using its SMILES string
        compounds = pcp.get_compounds(smiles, 'smiles')
        
        if compounds:
            # Get the first compound from the search results
            compound = compounds[0]
            
            # Find a common name from the list of synonyms. We'll look for "Tryptophan hydroxamate".
            product_name = "Tryptophan hydroxamate"
            is_common_name_found = False
            for synonym in compound.synonyms:
                if product_name.lower() == synonym.lower():
                    is_common_name_found = True
                    break
            
            # If the common name isn't found, fall back to the IUPAC name.
            if not is_common_name_found:
                 product_name = compound.iupac_name

            print(f"The name of the final product is: {product_name}")

        else:
            # Fallback name if the query fails
            product_name = "Tryptophan hydroxamate"
            print(f"Could not retrieve name from PubChem. Based on chemical analysis, the product is: {product_name}")

    except ImportError:
        print("Please install the 'pubchempy' library by running: pip install pubchempy")
    except Exception as e:
        print(f"An error occurred: {e}")
        print("Based on chemical analysis, the product is: Tryptophan hydroxamate")

find_chemical_name()