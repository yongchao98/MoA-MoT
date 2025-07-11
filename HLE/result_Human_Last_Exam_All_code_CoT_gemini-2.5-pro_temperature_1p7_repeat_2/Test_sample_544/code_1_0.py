# First, you may need to install the pubchempy library:
# pip install pubchempy

import pubchempy as pcp

def get_product_iupac_name():
    """
    This function determines the IUPAC name of the product from the reaction of
    methyl phenyl sulfoxide with triflic anhydride and trimethylsilyl cyanide.
    
    The reaction is a Pummerer rearrangement.
    
    Reaction breakdown:
    1. Reactants:
       - Methyl phenyl sulfoxide (C6H5-S(=O)-CH3)
       - Triflic anhydride ((CF3SO2)2O)
       - Trimethylsilyl cyanide ((CH3)3SiCN)
    2. Stoichiometry: The problem specifies 1 equivalent of each reagent.
       1 * C6H5-S(=O)-CH3 + 1 * (CF3SO2)2O + 1 * (CH3)3SiCN -> Product
    3. Mechanism:
       - The sulfoxide is activated by the triflic anhydride.
       - A proton is removed from the methyl group to form a sulfur ylide.
       - The ylide rearranges to a thionium ion ([C6H5-S=CH2]+).
       - The cyanide ion (CN-) from trimethylsilyl cyanide attacks the thionium ion.
    4. Product Structure: The final product is 2-(phenylthio)acetonitrile (C6H5-S-CH2-CN).
    5. Find IUPAC name: We use the product's SMILES string to find its name in the PubChem database.
    """
    
    # SMILES string for the product, 2-(phenylthio)acetonitrile
    product_smiles = "N#CCSc1ccccc1"
    
    print("Determining the IUPAC name for the product...")
    print(f"Product SMILES string: {product_smiles}")
    
    try:
        # Query PubChem for compounds matching the SMILES string
        compounds = pcp.get_compounds(product_smiles, 'smiles')
        
        if not compounds:
            # Fallback name if PubChem search fails
            product_name = "2-(phenylthio)acetonitrile"
            print("\nCould not find the compound in the PubChem database.")
            print(f"Based on chemical nomenclature rules, a correct name is: {product_name}")
            return
            
        # Get the first result
        product = compounds[0]
        
        # Get the official IUPAC name
        iupac_name = product.iupac_name
        
        if iupac_name:
            print(f"\nThe IUPAC name of the product is: {iupac_name}")
        else:
            # Fallback if the IUPAC name field is empty
            product_name = "2-(phenylthio)acetonitrile"
            print("\nIUPAC name not found in the database entry.")
            print(f"Based on chemical nomenclature rules, a correct name is: {product_name}")

    except Exception as e:
        print(f"\nAn error occurred while contacting the PubChem database: {e}")
        product_name = "2-(phenylthio)acetonitrile"
        print(f"Based on chemical nomenclature rules, a correct name is: {product_name}")

if __name__ == "__main__":
    get_product_iupac_name()