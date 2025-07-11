def get_reaction_product():
    """
    This function determines and prints the product of the specified chemical reaction.
    
    The reaction is a Directed ortho-Metalation followed by electrophilic quenching.
    1. N,N-diethyl-3-dimethylaminobenzamide has two directing groups (-CONEt2 and -NMe2).
    2. Both groups direct the base (sec-BuLi) to the C2 position.
    3. The resulting aryllithium is quenched with the electrophile (methyl iodide).
    4. A methyl group is added at the C2 position.
    """
    
    # Define the name of the final product based on the chemical analysis
    product_name = "N,N-diethyl-2-methyl-3-dimethylaminobenzamide"
    
    # Define the molecular formula of the product
    # Starting Material: C13H20N2O
    # Reagents add a methyl group (CH3) and remove a hydrogen (H). Net addition: CH2.
    # Final Formula: C(13+1)H(20+2)N2O = C14H22N2O
    product_formula = {
        'Carbon': 14,
        'Hydrogen': 22,
        'Nitrogen': 2,
        'Oxygen': 1
    }
    
    # Print the name of the compound
    print(f"The final compound obtained is: {product_name}")
    
    # Print the molecular formula information
    print(f"\nThe molecular formula is C{product_formula['Carbon']}H{product_formula['Hydrogen']}N{product_formula['Nitrogen']}O{product_formula['Oxygen']}.")
    
    print("\nThe count for each element in the final molecular formula is:")
    # Print each number in the "final equation" (elemental count)
    print(f"Carbon (C): {product_formula['Carbon']}")
    print(f"Hydrogen (H): {product_formula['Hydrogen']}")
    print(f"Nitrogen (N): {product_formula['Nitrogen']}")
    print(f"Oxygen (O): {product_formula['Oxygen']}")

# Execute the function to display the results
get_reaction_product()