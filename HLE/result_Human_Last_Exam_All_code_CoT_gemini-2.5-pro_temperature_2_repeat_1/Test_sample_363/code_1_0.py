def get_product_name():
    """
    This function determines and prints the IUPAC name of the reaction product.
    The reaction is an aza-Claisen rearrangement, a type of [3,3] sigmatropic shift.
    """
    
    # The chiral auxiliary on the nitrogen atom remains, but is now on a secondary amide (N-H).
    # Its structure includes the stereocenter (S) and locant 1.
    n_substituent = "N-((S)-1-phenylethyl)"

    # The rearrangement creates a new stereocenter at the alpha-carbon (position 2) of the propanamide chain.
    # The (S)-auxiliary directs the formation of the (R) stereocenter at position 2.
    new_stereo_and_locant = "(2R)-2-"
    
    # The migrating group, which becomes the substituent at position 2, is the original N-allyl group.
    # Its name, including its stereocenter (S) and locants (5, 1, 1), is maintained.
    migrating_substituent = "(((S)-5-methylcyclopent-1-en-1-yl)methyl)"
    
    # The parent acyl chain is a propanamide.
    parent_amide = "propanamide"
    
    # Assemble the full IUPAC name
    full_product_name = n_substituent + new_stereo_and_locant + migrating_substituent + parent_amide
    
    print("The IUPAC name of the final product is:")
    print(full_product_name)

get_product_name()