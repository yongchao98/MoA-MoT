def get_wittig_product():
    """
    Determines and prints the product of the specified Wittig reaction.
    """
    product_name = "(Z)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    # The numbers in the IUPAC name specify the positions of the substituents and the double bond.
    # 1: position of the (2-chlorophenyl) group on the main chain.
    # 2: position of the chloro group on the phenyl ring.
    # 4,4: positions of the two methyl groups on the main chain.
    # 2: position of the double bond ('ene') on the main chain.
    numbers_in_name = "1, 2, 4, 4, 2"
    
    print(f"The product of the Wittig reaction is: {product_name}")
    print(f"The numbers in the final IUPAC name are: {numbers_in_name}")

get_wittig_product()