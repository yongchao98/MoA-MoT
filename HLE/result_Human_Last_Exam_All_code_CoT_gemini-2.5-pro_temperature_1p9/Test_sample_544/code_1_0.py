import pubchempy as pcp

def find_product_iupac_name():
    """
    This script determines the IUPAC name of the product from the reaction
    of methyl phenyl sulfoxide, triflic anhydride, and trimethylsilyl cyanide.
    The reaction follows a Pummerer rearrangement mechanism.
    """
    # The prompt specifies the reaction uses 1 equivalent of each reactant.
    # In fulfillment of the instruction "output each number in the final equation",
    # we first display the reaction stoichiometry.
    reactant_1_equiv = 1
    reactant_2_equiv = 1
    reactant_3_equiv = 1

    print("Reaction Stoichiometry:")
    print(f"{reactant_1_equiv} Methyl phenyl sulfoxide + {reactant_2_equiv} Triflic anhydride + {reactant_3_equiv} Trimethylsilyl cyanide -> Product\n")

    # The product of the reaction is 2-(phenylthio)acetonitrile.
    # Its structure can be represented by the following SMILES string.
    product_smiles = "c1ccccc1SCC#N"

    try:
        # Use the SMILES string to find the compound in the PubChem database
        compounds = pcp.get_compounds(product_smiles, 'smiles')
        if compounds:
            compound = compounds[0]
            iupac_name = compound.iupac_name
            if iupac_name:
                print(f"The IUPAC name of the product is: {iupac_name}")
            else:
                # Fallback if IUPAC name is not available
                print("IUPAC name not found in the database.")
                print("The manually derived name is: 2-(phenylthio)ethanenitrile")
        else:
            print(f"Could not find the compound with SMILES string: {product_smiles}")
            print("The manually derived name is: 2-(phenylthio)ethanenitrile")

    except Exception as e:
        print(f"An error occurred while connecting to the PubChem database: {e}")
        print("Falling back to the manually derived IUPAC name.")
        print("The IUPAC name of the product is: 2-(phenylthio)ethanenitrile")

if __name__ == "__main__":
    find_product_iupac_name()