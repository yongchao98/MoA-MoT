def get_product_iupac_name():
    """
    This function returns the IUPAC name of the product from the described reaction.
    The reaction is a Pummerer reaction of methyl phenyl sulfoxide with triflic
    anhydride and trimethylsilyl cyanide.
    The final product is Ph-S-CH2-CN.
    """
    # The parent chain is acetonitrile (or ethanenitrile).
    # The substituent at position 2 is a phenylsulfanyl group.
    # The number '2' is part of the IUPAC name.
    product_name = "2-(Phenylsulfanyl)acetonitrile"
    print(f"The IUPAC name of the final product is: {product_name}")

get_product_iupac_name()