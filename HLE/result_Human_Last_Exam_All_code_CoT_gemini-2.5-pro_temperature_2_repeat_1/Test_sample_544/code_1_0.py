# The following script determines the IUPAC name of the reaction product.
# For this script to run, you first need to install the `pubchempy` library:
# pip install pubchempy

import pubchempy as pcp
import re
import sys

def find_product_iupac_name():
    """
    This function identifies the product of the specified reaction,
    retrieves its IUPAC name from PubChem, and prints the result.
    """
    # The reaction product is phenylthioacetonitrile.
    # Its structure can be represented by the SMILES string 'c1ccccc1SCC#N'.
    product_smiles = "c1ccccc1SCC#N"

    print(f"Querying database for product with SMILES: {product_smiles}")
    
    try:
        # Use the SMILES string to get compound information from PubChem.
        compounds = pcp.get_compounds(product_smiles, 'smiles')

        if not compounds:
            print("\nCould not find the compound in the PubChem database.")
            return

        # Get the first compound from the results.
        product_compound = compounds[0]

        # Get the IUPAC name.
        iupac_name = product_compound.iupac_name

        if iupac_name:
            print("\n--- Product Information ---")
            print(f"Final Product IUPAC Name: {iupac_name}")

            # Extract and print any numbers from the IUPAC name, as requested.
            # In this case, it corresponds to the locant number in the chemical name.
            numbers = re.findall(r'\d+', iupac_name)
            if numbers:
                print("\nNumber(s) found in the name:")
                for number in numbers:
                    print(number)
            print("-------------------------")

        else:
            print("IUPAC name could not be retrieved.")

    except ImportError:
        print("\nError: The 'pubchempy' library is not installed.", file=sys.stderr)
        print("Please install it using: pip install pubchempy", file=sys.stderr)
    except Exception as e:
        print(f"\nAn error occurred: {e}", file=sys.stderr)
        print("Please check your internet connection.", file=sys.stderr)

if __name__ == '__main__':
    find_product_iupac_name()