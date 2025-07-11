def get_iupac_name():
    """
    This function determines and prints the IUPAC name of the reaction product.
    
    Reactants:
    1. Methyl phenyl sulfoxide (C6H5-S(=O)-CH3)
    2. Triflic anhydride ((CF3SO2)2O) - Activator
    3. Trimethylsilyl cyanide ((CH3)3SiCN) - Nucleophile source
    
    Reaction: Pummerer reaction
    
    Mechanism:
    1. Sulfoxide is activated by triflic anhydride to form [C6H5-S(OTf)-CH3]+.
    2. A thionium ion intermediate, [C6H5-S=CH2]+, is formed.
    3. Cyanide from TMSCN attacks the CH2 carbon of the thionium ion.
    
    Product Structure: C6H5-S-CH2-CN
    
    IUPAC Naming:
    - Principal functional group: Nitrile (-CN) -> suffix "nitrile"
    - Parent alkane chain with the nitrile: CH2-CN -> Acetonitrile
    - Substituent on the acetonitrile: C6H5-S- group (phenylsulfanyl)
    - Position of the substituent: on the C2 atom of acetonitrile (C1 is the nitrile carbon)
    - Final Name: 2-(Phenylsulfanyl)acetonitrile
    """
    
    product_name = "2-(Phenylsulfanyl)acetonitrile"
    print(product_name)

get_iupac_name()