def get_product_iupac_name():
    """
    This function determines and prints the IUPAC name of the reaction product.
    The reaction is the complete substitution of halogens in 1,3-dibromo-2-iodobenzene
    with phenyl groups from excess phenyl magnesium bromide.
    """
    # The numbers for the positions of the substituents on the parent benzene ring.
    locants = "1,2,3"
    
    # The prefix 'tri' indicates three identical substituents.
    # The substituent is a 'phenyl' group.
    substituent_part = "triphenyl"
    
    # The parent molecule is benzene.
    parent_chain = "benzene"
    
    # Construct the final IUPAC name by combining the parts.
    # The format is locants-substituent_parent
    final_iupac_name = f"{locants}-{substituent_part}{parent_chain}"
    
    print("The IUPAC name of the final product is:")
    print(final_iupac_name)

# Execute the function to find the answer.
get_product_iupac_name()