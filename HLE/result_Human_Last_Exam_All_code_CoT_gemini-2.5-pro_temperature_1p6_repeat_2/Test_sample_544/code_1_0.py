def get_iupac_name():
    """
    This function determines and prints the IUPAC name of the reaction product.

    The reaction is a Pummerer rearrangement.
    1. Reactants: methyl phenyl sulfoxide (C6H5-S(=O)-CH3), triflic anhydride ((CF3SO2)2O),
       and trimethylsilyl cyanide ((CH3)3SiCN).
    2. The sulfoxide is activated by triflic anhydride, and after proton loss, it forms
       a thionium ion intermediate: [C6H5-S=CH2]+.
    3. The cyanide nucleophile (from TMSCN) attacks the electrophilic carbon of the
       thionium ion.
    4. The final product structure is C6H5-S-CH2-CN.
    5. Naming the product:
       - The principal functional group is the nitrile (-CN), which has higher priority
         than the sulfide (-S-).
       - The parent alkane is ethanenitrile (for the CH3-CN skeleton).
       - The nitrile carbon is C1, which means the substituent is on C2.
       - The substituent C6H5-S- is named 'phenylsulfanyl'.
       - Combining these parts gives the final IUPAC name.
    """
    iupac_name = "2-(phenylsulfanyl)ethanenitrile"
    print("The IUPAC name of the product is:")
    print(iupac_name)

get_iupac_name()