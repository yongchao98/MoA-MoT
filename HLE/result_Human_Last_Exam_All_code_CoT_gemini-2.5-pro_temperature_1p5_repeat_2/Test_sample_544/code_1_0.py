def get_iupac_name():
    """
    This function determines and prints the IUPAC name of the reaction product.
    
    The reaction is a Pummerer-type reaction:
    Reactants: methyl phenyl sulfoxide, triflic anhydride, trimethylsilyl cyanide
    Product: Ph-S-CH2-CN

    Naming follows IUPAC rules:
    1. Principal group: Nitrile (-CN)
    2. Parent chain: ethanenitrile (-CH2-CN)
    3. Substituent: phenylsulfanyl (Ph-S-) at position 2
    """
    
    # The systematic IUPAC name for the product Ph-S-CH2-CN
    iupac_name = "2-(Phenylsulfanyl)ethanenitrile"
    
    print(iupac_name)

if __name__ == "__main__":
    get_iupac_name()