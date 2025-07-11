def get_product_name():
    """
    This function returns the IUPAC name of the reaction product.
    The reaction is a stereoselective aza-Claisen rearrangement.
    """
    # The starting material is N-(((S)-5-methylcyclopent-1-en-1-yl)methyl)-N-((S)-1-phenylethyl)propionamide.
    # Step 1: LiHMDS forms the (Z)-enolate.
    # Step 2: Heating at 100°C causes a [3,3]-sigmatropic rearrangement.
    # The final product, after workup, is a carboxylic acid.
    
    # IUPAC Name determination:
    # Parent chain: propanoic acid
    # Stereocenter on the acid chain (C2): (S)
    # Substituent: A cyclopentyl ring attached to C2 of the acid.
    # The ring has a methylidene (=CH2) group and a methyl group.
    # Numbering the ring: C1 is the point of attachment, C2 has the methylidene, C4 has the methyl.
    # Stereocenters on the ring: C1 is (R), C4 is (S).
    # Full substituent name: ((1R,4S)-4-methyl-2-methylidenecyclopentyl)
    # Final product name combines these parts.
    
    product_name = "(2S)-2-((1R,4S)-4-methyl-2-methylidenecyclopentyl)propanoic acid"
    
    # Printing the numbers within the final name as requested by the prompt.
    # The numbers are: 2, 1, 4, 4, 2
    
    print(product_name)

get_product_name()
