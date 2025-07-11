def get_product_name():
    """
    This function returns the IUPAC name of the product from the reaction scheme.
    
    The reaction proceeds via two main stages:
    1. Knoevenagel condensation followed by dehydration to form diethyl (1-(ethoxycarbonyl)vinyl)phosphonate.
    2. A Michael-addition-Induced Ring Closure (MIRC) reaction. A thiolate adds to the vinylphosphonate,
       followed by an intramolecular aldol condensation and subsequent dehydration to yield a substituted
       dihydrothiophene ring.
    """
    
    # Structure analysis leads to a 4,5-dihydrothiophene core.
    # The principal functional group is an ethyl ester at position 4.
    # A diethoxyphosphoryl group is also at position 4.
    
    principal_group_locant = 4
    substituent_locant = 4
    
    # Based on the final dehydrated structure from the Michael-Aldol-Dehydration cascade.
    product_name = "ethyl {0}-(diethoxyphosphoryl)-4,5-dihydrothiophene-{1}-carboxylate".format(substituent_locant, principal_group_locant)
    
    print("The IUPAC name of the product is:")
    print(product_name)

get_product_name()