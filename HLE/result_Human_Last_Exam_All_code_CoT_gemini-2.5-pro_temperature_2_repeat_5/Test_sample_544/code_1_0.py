def get_pummerer_product_name():
    """
    This function determines the IUPAC name of the product from the specified reaction.
    
    The reaction is a Pummerer reaction variant involving:
    - Methyl phenyl sulfoxide
    - Triflic anhydride (activator)
    - Trimethylsilyl cyanide (nucleophile)

    The reaction proceeds through a thionium ion intermediate ([Ph-S=CH2]+),
    which is then attacked by the cyanide nucleophile.

    The final product has the structure Ph-S-CH2-CN.

    Its IUPAC name is derived by treating 'acetonitrile' as the parent
    compound and 'phenylthio' as the substituent.
    """
    
    # Define the components of the final product name
    substituent = "Phenylthio"
    parent_chain = "acetonitrile"
    position = 2 # The substituent is on carbon-2 of the acetonitrile chain
    
    # Construct the IUPAC name
    iupac_name = f"{position}-({substituent}){parent_chain}"
    
    # The name is 2-(Phenylthio)acetonitrile
    print(iupac_name)

get_pummerer_product_name()