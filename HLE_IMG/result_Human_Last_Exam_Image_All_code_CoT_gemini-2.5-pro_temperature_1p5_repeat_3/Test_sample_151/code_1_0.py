# First, ensure you have the pubchempy library installed:
# pip install pubchempy

import pubchempy as pcp

def get_iupac_name_from_smiles(smiles_string):
    """
    Finds a compound on PubChem using its SMILES string and prints its IUPAC name.
    """
    try:
        # Search for compounds with the given SMILES string
        compounds = pcp.get_compounds(smiles_string, 'smiles')
        if compounds:
            # Take the first result
            compound = compounds[0]
            print(f"SMILES: {smiles_string}")
            # The 'best' name is often the IUPAC name, but let's be explicit
            iupac_name = compound.iupac_name
            if iupac_name:
                print(f"IUPAC Name: {iupac_name}")
            else:
                print("IUPAC name not found, but other names exist.")
        else:
            print("No compound found for the given SMILES string.")
    except Exception as e:
        print(f"An error occurred: {e}")

# The SMILES string for the final product, ethyl 2,5-dihydrothiophene-3-carboxylate
product_smiles = "S1CC(C(=O)OCC)=CC1"

get_iupac_name_from_smiles(product_smiles)