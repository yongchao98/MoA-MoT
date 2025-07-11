def get_iupac_name():
    """
    This function returns the IUPAC name of the product from the reaction of
    methyl phenyl sulfoxide with triflic anhydride and trimethylsilyl cyanide.
    
    The reaction is a Pummerer rearrangement.
    1. Reactants: Methyl phenyl sulfoxide (Ph-S(=O)-Me), Triflic anhydride (Tf2O), Trimethylsilyl cyanide (TMSCN).
    2. Activation: Sulfoxide attacks Tf2O to form [Ph-S(OTf)-Me]+.
    3. Thionium ion formation: Deprotonation of the methyl group yields the electrophilic thionium ion, [Ph-S=CH2]+.
    4. Nucleophilic attack: Cyanide (from TMSCN) attacks the CH2 carbon.
    5. Product: The final product is Ph-S-CH2-CN.
    
    IUPAC Naming of Ph-S-CH2-CN:
    - Principal functional group: Nitrile (-CN)
    - Parent chain: Ethanenitrile (or acetonitrile)
    - Substituent: The C6H5-S- group is called 'phenylthio'.
    - Position: It is on the second carbon (C2) of the ethanenitrile chain.
    - Final Name: 2-(phenylthio)acetonitrile
    """
    product_name = "2-(phenylthio)acetonitrile"
    print(product_name)

get_iupac_name()