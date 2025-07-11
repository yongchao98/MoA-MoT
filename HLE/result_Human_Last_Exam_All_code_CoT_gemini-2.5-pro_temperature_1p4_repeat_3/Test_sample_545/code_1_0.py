def get_iupac_name():
    """
    This function returns the IUPAC name of the major product.
    The reaction is a sulfoxide syn-elimination.
    Reactant: ((2-((2-methylbut-3-en-2-yl)oxy)ethyl)sulfinyl)benzene
    Conditions: Heat at 180 C in decalin with NaHCO3
    Product Structure: CH2=CH-O-C(CH3)2-CH=CH2
    IUPAC Name is determined using replacement nomenclature.
    """
    # The numbers in the name are locants indicating positions of substituents and functional groups.
    # 4,4-dimethyl: Two methyl groups at position 4
    # 3-oxa: An oxygen atom replaces the carbon at position 3
    # hexa-1,5-diene: A 6-carbon chain with double bonds at positions 1 and 5
    product_name = "4,4-dimethyl-3-oxahexa-1,5-diene"
    
    # Printing the final name, which includes all the required numbers.
    print(product_name)

get_iupac_name()