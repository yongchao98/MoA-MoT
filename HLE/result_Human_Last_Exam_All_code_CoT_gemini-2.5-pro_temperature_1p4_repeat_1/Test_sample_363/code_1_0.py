def get_product_name():
    """
    This function returns the IUPAC name of the product from the described reaction.
    The reaction is a Myers asymmetric aza-Claisen rearrangement.
    """
    
    # The reaction yields a single major diastereomer.
    # The naming components are determined by the stereoselective transformation.
    
    # Stereochemistry at the alpha-carbon of the acyl chain, controlled by the auxiliary.
    acyl_chain_stereo = "(2R)"
    
    # The chiral auxiliary on the nitrogen atom.
    n_substituent_auxiliary = "N-((S)-1-phenylethyl)"
    
    # The location of the main substituent on the acyl chain.
    acyl_chain_locant = "2"
    
    # The name and stereochemistry of the rearranged substituent.
    # Stereochemistry is determined by the reaction mechanism and existing chirality.
    rearranged_substituent = "((1R,4S)-4-methyl-2-methylidenecyclopentyl)"
    
    # The base name of the new acyl group.
    acyl_chain_base = "propanamide"
    
    # Assemble the full IUPAC name.
    product_name = f"{acyl_chain_stereo}-{n_substituent_auxiliary}-{acyl_chain_locant}-{rearranged_substituent}{acyl_chain_base}"
    
    print(product_name)

get_product_name()