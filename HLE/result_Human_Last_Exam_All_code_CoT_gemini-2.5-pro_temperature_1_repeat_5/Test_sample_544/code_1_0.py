def get_iupac_name():
    """
    This function determines and prints the IUPAC name of the product from the described reaction.
    """
    # The reaction is a Pummerer reaction.
    # Reactants: Methyl phenyl sulfoxide, Triflic anhydride, Trimethylsilyl cyanide.
    # Product structure: C6H5-S-CH2-CN

    # IUPAC naming convention:
    # Principal functional group: Nitrile (-CN)
    # Parent name: Acetonitrile (common) or Ethanenitrile (systematic)
    # Substituent: Phenylthio (C6H5-S-)
    # Position of substituent: 2

    # The number in the IUPAC name specifies the location of the substituent on the main carbon chain.
    locant = 2
    substituent = "phenylthio"
    parent_chain = "acetonitrile"

    final_name = f"{locant}-({substituent}){parent_chain}"

    print(f"The IUPAC name of the product is: {final_name}")
    print(f"The number in the final name is: {locant}")

get_iupac_name()