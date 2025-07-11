def get_iupac_name_of_product():
    """
    This function determines and prints the IUPAC name for the product of a specific
    chemical reaction.

    The reaction involves:
    - Reactant 1: methyl phenyl sulfoxide
    - Reactant 2: 1 equivalent of triflic anhydride (activator)
    - Reactant 3: 1 equivalent of trimethylsilyl cyanide (nucleophile)

    This is a Pummerer-type reaction, which transforms a sulfoxide into an
    alpha-functionalized sulfide.

    The resulting product is 2-(phenylthio)ethanenitrile.
    The number '2' in the name indicates the position of the 'phenylthio'
    substituent on the 'ethanenitrile' parent chain.
    """

    # The IUPAC name of the final product, determined through chemical analysis.
    product_iupac_name = "2-(phenylthio)ethanenitrile"

    print("The IUPAC name of the reaction product is:")
    print(product_iupac_name)

# Execute the function to print the result.
get_iupac_name_of_product()