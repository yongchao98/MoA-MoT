def get_product_name():
    """
    This function constructs and prints the IUPAC name of the reaction product.
    The reaction is a stereoselective aza-Claisen rearrangement.
    """
    
    # Stereochemistry of the main chain's new chiral center (C2)
    # The (S)-phenylethyl auxiliary directs the formation of the (R) stereocenter.
    c2_config = "(2R)"
    
    # Stereochemistry of the substituent ring
    # The new chiral center (C2) and the existing one (C4) are determined by the reaction mechanism.
    ring_config = "(2S,4S)"
    
    # Substituent group name derived from the rearranged allyl moiety
    substituent_group = f"-{ring_config}-4-methyl-1-methylidenecyclopentan-2-yl"
    
    # Chiral auxiliary on the Nitrogen atom
    n_substituent = "N-((S)-1-phenylethyl)"
    
    # Main chain name
    parent_chain = "propanamide"
    
    # Assembling the final IUPAC name
    # The structure is a propanamide with two substituents at position 2.
    # The numbers in the name refer to the carbon positions as per IUPAC rules.
    final_name = (
        f"{c2_config}-2-methyl-2-({substituent_group.strip('-')})"
        f"-{n_substituent}{parent_chain}"
    )

    print("The IUPAC name of the product is:")
    print(final_name)
    
    # Output each number in the final equation (name)
    print("\nBreaking down the numbers in the name:")
    print("2R: Stereoconfiguration at carbon 2 of the propanamide chain.")
    print("2: Position of the methyl group on the propanamide chain.")
    print("2: Position of the cyclopentyl substituent on the propanamide chain.")
    print("2S: Stereoconfiguration at carbon 2 of the cyclopentyl ring.")
    print("4S: Stereoconfiguration at carbon 4 of the cyclopentyl ring.")
    print("4: Position of the methyl group on the cyclopentyl ring.")
    print("1: Position of the methylidene group on the cyclopentyl ring.")
    print("2: Attachment point of the cyclopentyl ring to the main chain.")
    print("S: Stereoconfiguration of the 1-phenylethyl group.")
    print("1: Position of the phenyl group on the ethyl chain of the auxiliary.")


get_product_name()
