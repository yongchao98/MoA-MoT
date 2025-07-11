def get_iupac_name():
    """
    This function returns the IUPAC name of the major product from the described reaction.
    
    The reaction is a sulfoxide pyrolysis (syn-elimination).
    Reactant: ((2-((2-methylbut-3-en-2-yl)oxy)ethyl)sulfinyl)benzene
    Structure: Ph-S(=O)-CH2-CH2-O-C(CH3)(CH=CH2)-CH3
    
    The elimination reaction breaks the C-S bond and a beta-C-H bond,
    forming an alkene and phenylsulfenic acid.
    
    Product Structure: CH2=CH-O-C(CH3)(CH=CH2)-CH3
    
    IUPAC Naming of the Product:
    1. Parent chain: but-1-ene (C1H2=C2H-C3H-C4H3)
    2. Substituents on the parent chain are at position 3.
    3. The substituents are a methyl group (-CH3) and an ethenyloxy group (-O-CH=CH2).
    4. Alphabetizing the substituents gives the name.
    """
    
    # The IUPAC name is composed of several parts with specific numbering.
    substituent_1 = "ethenyloxy"
    position_1 = 3
    
    substituent_2 = "methyl"
    position_2 = 3
    
    parent_chain = "but"
    double_bond_position = 1
    parent_suffix = "ene"
    
    # Construct the final name string, ensuring all numbers are included.
    # The name is 3-(ethenyloxy)-3-methylbut-1-ene
    final_name = (
        f"{position_1}-({substituent_1})-"
        f"{position_2}-{substituent_2}{parent_chain}-"
        f"{double_bond_position}-{parent_suffix}"
    )
    
    print(final_name)

get_iupac_name()