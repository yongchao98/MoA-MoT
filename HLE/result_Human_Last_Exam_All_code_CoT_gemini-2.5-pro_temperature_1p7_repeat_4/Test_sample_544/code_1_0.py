def get_reaction_product_name():
    """
    This function determines the IUPAC name of the product from the reaction of
    methyl phenyl sulfoxide with triflic anhydride and trimethylsilyl cyanide.

    The reaction proceeds via a Pummerer-type mechanism:
    1.  The sulfoxide is activated by the electrophilic triflic anhydride.
    2.  A thionium ion intermediate ([Ph-S=CH2]+) is formed via elimination.
    3.  The cyanide nucleophile (from TMSCN) attacks the thionium ion.

    The final product has the structure C6H5-S-CH2-CN.

    The IUPAC name is determined as follows:
    - Principal functional group: -CN (nitrile)
    - Parent chain: Acetonitrile (or ethanenitrile)
    - Substituent: A 'phenylsulfanyl' group (C6H5-S-) is on the carbon
      at position 2.

    The resulting name is 2-(phenylsulfanyl)acetonitrile.
    """
    iupac_name = "2-(phenylsulfanyl)acetonitrile"
    print(f"The IUPAC name of the product is: {iupac_name}")

if __name__ == '__main__':
    get_reaction_product_name()