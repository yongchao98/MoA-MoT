# First, ensure you have the pubchempy library installed:
# pip install pubchempy

import pubchempy as pcp

def find_iupac_name_for_pummerer_product():
    """
    This script determines the IUPAC name of the product from the reaction of
    methyl phenyl sulfoxide with triflic anhydride and trimethylsilyl cyanide.

    1. The reaction is a Pummerer rearrangement.
    2. The sulfoxide group is activated and replaced, and the cyanide nucleophile
       attacks the alpha-carbon.
    3. The resulting product is (phenylthio)acetonitrile, which can be
       represented by the SMILES string 'N#CCSc1ccccc1'.
    4. This script uses the SMILES string to find the compound's official
       IUPAC name from the PubChem database.
    """
    try:
        # The SMILES string for the product, 2-(phenylthio)ethanenitrile.
        product_smiles = "N#CCSc1ccccc1"

        # Use PubChemPy to get the compound from the SMILES string.
        compounds = pcp.get_compounds(product_smiles, 'smiles')

        if compounds:
            # Get the first compound found.
            compound = compounds[0]
            # Retrieve the official IUPAC name.
            iupac_name = compound.iupac_name

            if iupac_name:
                print("The IUPAC name of the product is:")
                # The IUPAC name is 2-(phenylsulfanyl)ethanenitrile.
                # The prompt asks to output each number in the final 'equation'/name.
                # We will print the full name string, which includes the number '2'.
                print(iupac_name)
            else:
                # Fallback in case the IUPAC name is not available.
                print("IUPAC name not found in the database. A valid name is 2-(phenylsulfanyl)ethanenitrile.")

        else:
            print(f"Could not find a compound with SMILES: {product_smiles}")
            print("The expected product is 2-(phenylsulfanyl)ethanenitrile.")

    except Exception as e:
        print(f"An error occurred while communicating with PubChem: {e}")
        print("Based on chemical principles, the IUPAC name is 2-(phenylsulfanyl)ethanenitrile.")

if __name__ == "__main__":
    find_iupac_name_for_pummerer_product()