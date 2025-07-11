import pubchempy as pcp

def get_product_iupac_name():
    """
    This function determines the IUPAC name of the product from the reaction of
    methyl phenyl sulfoxide, triflic anhydride, and trimethylsilyl cyanide.

    The reaction is a Pummerer rearrangement, where the alpha-methyl group
    of the sulfoxide is functionalized by a cyanide nucleophile.

    Reactant (SMILES): CS(=O)c1ccccc1
    Product structure: Ph-S-CH2-CN
    Product SMILES: c1ccccc1SCC#N

    This script retrieves the IUPAC name for the product from PubChem.
    """
    product_smiles = "c1ccccc1SCC#N"

    print(f"Querying PubChem for the compound with SMILES: {product_smiles}")

    try:
        # Search PubChem by the SMILES string
        compounds = pcp.get_compounds(product_smiles, 'smiles')
        if compounds:
            # Get the first result
            compound = compounds[0]
            iupac_name = compound.iupac_name
            # According to the prompt "output each number in the final equation!",
            # we will print the full name, including the locant '2'.
            print("\nThe IUPAC name of the product is:")
            print(iupac_name)
        else:
            # Fallback name if not found in the database
            print("Compound not found in PubChem. The deduced IUPAC name is:")
            print("2-(phenylsulfanyl)acetonitrile")

    except Exception as e:
        print(f"An error occurred while querying PubChem: {e}")
        print("The deduced IUPAC name is:")
        print("2-(phenylsulfanyl)acetonitrile")

if __name__ == "__main__":
    get_product_iupac_name()