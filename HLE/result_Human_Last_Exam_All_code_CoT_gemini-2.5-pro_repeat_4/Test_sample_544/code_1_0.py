def get_iupac_name():
    """
    This function returns the IUPAC name of the product from the described reaction.
    
    The reaction is the Pummerer rearrangement of methyl phenyl sulfoxide.
    Reactants:
    1. Methyl phenyl sulfoxide
    2. Triflic anhydride (activator)
    3. Trimethylsilyl cyanide (nucleophile source)
    
    The reaction proceeds via an activated sulfonium salt, which forms a thionium ion.
    The thionium ion is then trapped by the cyanide nucleophile.
    
    Product structure: Ph-S-CH2-CN
    
    IUPAC Naming:
    - Principal group: Nitrile (-CN) -> Parent: acetonitrile
    - Substituent: Phenylsulfanyl group (Ph-S-)
    - Position: Attached to carbon-2 of the acetonitrile backbone.
    """
    product_name = "2-(Phenylsulfanyl)acetonitrile"
    print(f"The IUPAC name of the product is: {product_name}")

get_iupac_name()