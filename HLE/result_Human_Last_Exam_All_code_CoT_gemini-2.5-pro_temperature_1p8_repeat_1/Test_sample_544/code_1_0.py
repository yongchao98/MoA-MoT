# To run this code, you may need to install the pubchempy library first:
# pip install pubchempy

import pubchempy as pcp

def find_reaction_product_name():
    """
    This function determines the product of the Pummerer reaction between
    methyl phenyl sulfoxide, triflic anhydride, and trimethylsilyl cyanide,
    and then retrieves its IUPAC name.
    """

    # Based on the Pummerer reaction mechanism, the product is an alpha-cyano sulfide.
    # Its structure is C6H5-S-CH2-CN.
    # We can represent this structure using a SMILES string.
    product_smiles = "N#CCSc1ccccc1"
    
    product_iupac_name = "Not found"

    try:
        # Use the SMILES string to find the compound in the PubChem database.
        compounds = pcp.get_compounds(product_smiles, 'smiles')
        if compounds:
            # Retrieve the standard IUPAC name from the database entry.
            product_iupac_name = compounds[0].iupac_name
        else:
            # Fallback name based on IUPAC nomenclature rules if not found.
            product_iupac_name = "2-(phenylsulfanyl)acetonitrile"
            
    except Exception as e:
        print(f"An error occurred (ensure pubchempy is installed and you have an internet connection): {e}")
        # Provide a fallback name in case of an error.
        product_iupac_name = "2-(phenylsulfanyl)acetonitrile"

    # Per the instructions, here is the final equation with the numbers (stoichiometric coefficients).
    print("Reaction Equation:")
    print("1 methyl phenyl sulfoxide + 1 triflic anhydride + 1 trimethylsilyl cyanide -> Product")
    print("\n")
    print("The IUPAC name of the main organic product is:")
    print(product_iupac_name)

# Execute the function to find and print the name.
find_reaction_product_name()