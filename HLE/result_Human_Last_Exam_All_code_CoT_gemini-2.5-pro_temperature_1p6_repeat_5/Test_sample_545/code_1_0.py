def get_iupac_name():
    """
    This function returns the IUPAC name of the major product from the described reaction.
    
    The reaction sequence is:
    1. Thermal sulfoxide elimination to form an allyl vinyl ether.
    2. [3,3]-Sigmatropic Claisen rearrangement of the ether to form a gamma,delta-unsaturated aldehyde.
    
    The final product is 3,3-dimethylpent-4-enal.
    The numbers in the name are locants:
    - 3,3: Indicate the position of the two methyl groups on the third carbon of the main chain.
    - 4: Indicates the starting position of the double bond on the fourth carbon of the main chain.
    """
    
    # The IUPAC name of the final product.
    final_product_name = "3,3-dimethylpent-4-enal"
    
    # Print the name as requested. The numbers 3, 3, and 4 are part of the name.
    print(final_product_name)

get_iupac_name()