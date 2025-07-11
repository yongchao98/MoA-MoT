# First, ensure you have the pubchempy library installed:
# pip install pubchempy

import sys
try:
    import pubchempy as pcp
except ImportError:
    print("The required library 'pubchempy' is not found.")
    print("Please install it by running: pip install pubchempy")
    sys.exit(1)

def find_iupac_name_from_smiles(smiles_string):
    """
    Finds the IUPAC name for a given SMILES string using the PubChem database.
    """
    try:
        # Search for compounds using the SMILES string
        compounds = pcp.get_compounds(smiles_string, 'smiles')
        
        if not compounds:
            return "Compound not found in PubChem database."
            
        # Take the first result
        compound = compounds[0]
        
        # Return the IUPAC name if available
        return compound.iupac_name if compound.iupac_name else "IUPAC name not available."
        
    except Exception as e:
        return f"An error occurred during the PubChem search: {e}"

if __name__ == "__main__":
    # The reaction of methyl phenyl sulfoxide with triflic anhydride and TMSCN
    # produces (phenylsulfanyl)acetonitrile.
    # Its SMILES string is 'c1ccccc1SCC#N'.
    product_smiles = 'c1ccccc1SCC#N'
    
    # Get the IUPAC name
    iupac_name = find_iupac_name_from_smiles(product_smiles)
    
    print(f"The IUPAC name of the product is:")
    print(iupac_name)
