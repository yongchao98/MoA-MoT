# First, you may need to install the pubchempy library.
# You can do this by running the following command in your terminal:
# pip install pubchempy

import pubchempy as pcp

def find_chemical_name(common_name):
    """
    Searches for a chemical by its common name using PubChem
    and prints its IUPAC name.
    """
    try:
        # Search for the compound by its common name
        compounds = pcp.get_compounds(common_name, 'name')
        if compounds:
            # Get the first search result
            product = compounds[0]
            # Retrieve the IUPAC name
            iupac_name = product.iupac_name
            print(f"The common name of the product is: {common_name}")
            print(f"The IUPAC name of the product is: {iupac_name}")
        else:
            print(f"Could not find a compound with the name: {common_name}")
    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an internet connection and the 'pubchempy' library is installed.")

# The product identified from the reaction scheme is Tryptophan hydroxamate.
# Let's use this name to find the formal IUPAC name.
product_common_name = "Tryptophan hydroxamate"
find_chemical_name(product_common_name)
