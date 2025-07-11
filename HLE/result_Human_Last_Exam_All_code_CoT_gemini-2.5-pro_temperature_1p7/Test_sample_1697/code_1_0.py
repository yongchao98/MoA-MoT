def get_reaction_product():
    """
    This function determines and prints the name of the product from the specified reaction.
    
    Reaction: N,N-diethyl-3-dimethylaminobenzamide reacts first with 
              sec-BuLi and TMEDA in THF and then with methyl iodide.
    
    Step 1: Directed ortho-Metalation.
    The N,N-diethylamide group at C-1 and the dimethylamino group at C-3 both direct
    the lithiation by sec-BuLi to the C-2 position.
    
    Step 2: Electrophilic Quench.
    The resulting aryllithium intermediate is a strong nucleophile that attacks
    the electrophile, methyl iodide (CH3-I).
    
    Result: A methyl group (-CH3) is added at the C-2 position.
    """
    
    starting_material = "N,N-diethyl-3-dimethylaminobenzamide"
    substituent_added = "2-methyl"
    
    # Constructing the final product name
    product_name = "N,N-diethyl-2-methyl-3-dimethylaminobenzamide"
    
    print(f"The starting material is: {starting_material}")
    print(f"The reaction is a directed ortho-metalation followed by methylation.")
    print(f"The final product obtained is: {product_name}")

get_reaction_product()