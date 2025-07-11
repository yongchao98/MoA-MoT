def get_reaction_product():
    """
    This function determines and prints the product of the specified multi-step organic reaction.
    The reaction is a Directed ortho-Metalation followed by an electrophilic quench.
    """

    # The starting material is N,N-diethyl-3-dimethylaminobenzamide.
    # The reaction sequence is:
    # 1) sec-BuLi, TMEDA, THF
    # 2) CH3I

    # The most acidic proton is at the C2 position, ortho to the very strong
    # amide directing group and also ortho to the dimethylamino directing group.
    # Lithiation occurs at C2.
    
    # Quenching with methyl iodide adds a methyl group at the C2 position.
    
    # The name of the final product is constructed accordingly.
    # The numbers in the name are the locants on the benzene ring (3 and 2)
    # and the designators for the substituents on the nitrogen atom (N,N).
    product_name = "N,N-diethyl-3-(dimethylamino)-2-methylbenzamide"
    
    print(f"The compound obtained from the reaction is: {product_name}")

# Execute the function to find the product
get_reaction_product()