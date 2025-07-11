def get_iupac_name():
    """
    This function determines and returns the IUPAC name of the reaction product.

    The reaction is a Pummerer-type rearrangement.
    Reactants:
    1. Methyl phenyl sulfoxide (C6H5-S(=O)-CH3)
    2. Triflic anhydride ((CF3SO2)2O)
    3. Trimethylsilyl cyanide ((CH3)3SiCN)

    Product Structure: C6H5-S-CH2-CN

    IUPAC Naming:
    - Principal functional group: Nitrile (-CN) -> suffix "ethanenitrile" for the 2-carbon chain.
    - Substituent at position 2: Phenylsulfanyl group (C6H5-S-).
    - The number '2' is the locant indicating the position of the substituent.
    """
    iupac_name = "2-(Phenylsulfanyl)ethanenitrile"
    print("The IUPAC name of the product is:")
    print(iupac_name)

if __name__ == '__main__':
    get_iupac_name()