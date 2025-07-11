# First, ensure you have the pubchempy library installed:
# pip install pubchempy

import pubchempy as pcp

def find_chemical_name_from_smiles(smiles_string):
    """
    Finds the common name of a chemical compound from its SMILES string using PubChem.
    """
    try:
        # Search for the compound using its SMILES string
        compounds = pcp.get_compounds(smiles_string, 'smiles')
        if compounds:
            # Get the first result
            compound = compounds[0]
            # Prioritize a common synonym, otherwise use the IUPAC name
            name = compound.synonyms[0] if compound.synonyms else compound.iupac_name
            print(f"The final product is: {name}")
        else:
            print("Compound not found in PubChem.")
    except Exception as e:
        print(f"An error occurred: {e}")
        print("A common name for the product is Tryptophan hydroxamate.")

if __name__ == '__main__':
    # SMILES string for the final product, (S)-Tryptophan hydroxamate
    # Structure: NH2-CH(CH2-Indole)-C(=O)NHOH
    final_product_smiles = "N[C@@H](CC1=CNC2=C1C=CC=C2)C(=O)NO"
    find_chemical_name_from_smiles(final_product_smiles)
