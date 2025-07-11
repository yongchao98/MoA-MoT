def get_product_iupac_name():
    """
    This script determines and prints the IUPAC name of the product resulting from the reaction
    of methyl phenyl sulfoxide with triflic anhydride and trimethylsilyl cyanide.
    
    The reaction is a Pummerer-type reaction.
    - Reactant: C6H5-S(=O)-CH3 (methyl phenyl sulfoxide)
    - Reagents: Triflic anhydride (activator) and TMSCN (nucleophile source)
    - Product Structure: C6H5-S-CH2-CN
    
    This script will now construct the IUPAC name for this product.
    """

    # According to IUPAC rules, the nitrile group (-CN) has priority over the sulfide (-S-).
    # Therefore, the parent compound is based on the nitrile.
    # The parent chain is CH3-CN, which is named acetonitrile.
    # In naming, the nitrile carbon is C1 and the adjacent carbon is C2.
    parent_chain = "acetonitrile"

    # The group C6H5-S- (phenyl group bonded to sulfur) is the substituent.
    # This substituent is named 'phenylthio'.
    substituent = "phenylthio"

    # The 'phenylthio' group is attached to the second carbon (C2) of the acetonitrile chain.
    # Therefore, the numerical locant for the substituent is 2.
    locant = 2

    # The final IUPAC name is assembled from these parts.
    # Format: locant-(substituent)parent_chain
    # Parentheses are used for the complex substituent 'phenylthio'.
    final_iupac_name = f"{locant}-({substituent}){parent_chain}"

    print("The IUPAC name for the product is constructed as follows:")
    print(f" - Parent Chain: {parent_chain}")
    print(f" - Substituent Group: {substituent}")
    print(f" - Position (Locant) of the Substituent: {locant}")
    print("\nFinal IUPAC Name:")
    print(final_iupac_name)

# Execute the function to find and print the name.
get_product_iupac_name()