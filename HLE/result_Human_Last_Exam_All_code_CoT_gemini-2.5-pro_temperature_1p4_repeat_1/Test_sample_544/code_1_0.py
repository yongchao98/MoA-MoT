import pubchempy as pcp
import sys

def find_iupac_name_of_reaction_product():
    """
    This script determines the IUPAC name of the product from the reaction of
    methyl phenyl sulfoxide with triflic anhydride and trimethylsilyl cyanide.

    The reaction is a Pummerer rearrangement, yielding a product whose structure
    can be represented by the SMILES string 'N#CCSc1ccccc1'.

    This script queries the PubChem database to find the IUPAC name for this structure.
    """
    # SMILES string for the product, which is 2-(phenylthio)acetonitrile.
    product_smiles = "N#CCSc1ccccc1"

    try:
        # Search for the compound in the PubChem database using its SMILES string.
        compounds = pcp.get_compounds(product_smiles, 'smiles')

        if not compounds:
            # Fallback name if PubChem is unavailable or the compound is not found.
            # This name is systematically derived.
            # -CN is a nitrile. The two-carbon chain is ethanenitrile.
            # The C6H5-S- group (phenylthio) is on carbon 2.
            iupac_name = "2-(phenylthio)ethanenitrile"
            print(f"Could not find the compound in PubChem. A valid systematic IUPAC name is: {iupac_name}")
            return

        # Extract the primary compound from the search results.
        product_compound = compounds[0]

        # Get the IUPAC name from the compound properties.
        iupac_name = product_compound.iupac_name

        print(f"The SMILES representation of the product is: {product_smiles}")
        print(f"The IUPAC name of the product is: {iupac_name}")

    except Exception as e:
        # Handle potential network errors or other issues with the library.
        print(f"An error occurred: {e}", file=sys.stderr)
        print("Please ensure you have a working internet connection and have installed the 'pubchempy' library (`pip install pubchempy`).")
        print("A valid IUPAC name for the product is '2-(phenylthio)acetonitrile'.")


if __name__ == "__main__":
    find_iupac_name_of_reaction_product()
