def get_iupac_name():
    """
    This function returns the IUPAC name of the product from the described reaction.
    
    The reaction is a Pummerer reaction.
    1. Reactants: methyl phenyl sulfoxide (Ph-S(=O)-Me), triflic anhydride ((Tf)2O), and trimethylsilyl cyanide (TMSCN).
    2. Activation: The sulfoxide oxygen attacks the triflic anhydride.
    3. Elimination: A proton is lost from the methyl group to form a thionium ion intermediate (Ph-S+=CH2).
    4. Nucleophilic attack: The cyanide ion (from TMSCN) attacks the CH2 carbon.
    5. Product: The resulting structure is Ph-S-CH2-CN.
    
    The IUPAC name is determined as follows:
    - Principal functional group: Nitrile (-CN).
    - Parent alkane: Ethane. So the parent name is ethanenitrile.
    - Substituent: The C6H5-S- group is called a "phenylthio" group.
    - Position: The phenylthio group is on carbon 2 of the ethanenitrile chain.
    """
    
    # Define the components of the IUPAC name
    position = 2
    substituent = "phenylthio"
    parent_chain = "ethanenitrile"
    
    # Construct the final name
    # The final name is 2-(phenylthio)ethanenitrile
    iupac_name = f"{position}-({substituent}){parent_chain}"
    
    print(iupac_name)

get_iupac_name()