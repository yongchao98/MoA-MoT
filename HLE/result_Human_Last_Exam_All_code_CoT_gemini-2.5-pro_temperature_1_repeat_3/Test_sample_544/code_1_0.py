# This script uses the pubchempy library to find the IUPAC name of a chemical.
# If you don't have it installed, you can run this command in your terminal:
# pip install pubchempy

import pubchempy as pcp

def find_product_iupac_name(common_name):
    """
    Searches PubChem for a compound by its common name and prints its IUPAC name.
    The reaction is the Pummerer reaction of methyl phenyl sulfoxide with
    triflic anhydride and trimethylsilyl cyanide.
    """
    print("Analyzing the Pummerer reaction...")
    print("---------------------------------")
    print("Reactant 1: Methyl phenyl sulfoxide")
    print("Reactant 2: Triflic anhydride (activator)")
    print("Reactant 3: Trimethylsilyl cyanide (nucleophile source)")
    print(f"\nThe predicted product has the structure C6H5-S-CH2-CN, commonly known as '{common_name}'.")
    print("Querying PubChem database for the official IUPAC name...")
    
    try:
        # Get a list of compounds matching the name
        compounds = pcp.get_compounds(common_name, 'name')
        
        if not compounds:
            # Fallback to chemically derived name if search fails
            derived_name = "2-(phenylthio)ethanenitrile"
            print(f"\nCould not find compound '{common_name}' in PubChem.")
            print(f"A chemically correct IUPAC name is: {derived_name}")
            return

        # Use the first and most relevant result
        product = compounds[0]
        
        # Retrieve the official IUPAC name
        iupac_name = product.iupac_name
        
        print("\n--- Search Result ---")
        print(f"The IUPAC name of the product is: {iupac_name}")

    except Exception as e:
        derived_name = "2-(phenylthio)ethanenitrile"
        print(f"\nAn error occurred during the search: {e}")
        print(f"A chemically correct IUPAC name is: {derived_name}")

# The common name for the product C6H5-S-CH2-CN
product_common_name = "(phenylthio)acetonitrile"
find_product_iupac_name(product_common_name)