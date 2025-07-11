def chemical_reaction_synthesis():
    """
    This script determines the product of a two-step organic synthesis reaction.
    It identifies the reaction type, analyzes the regioselectivity, and names the final product.
    """
    
    # Define the reaction components
    starting_material = "N,N-diethyl-3-dimethylaminobenzamide"
    reagents_step1 = "sec-BuLi, TMEDA, THF"
    reagents_step2 = "methyl iodide (CH3I)"
    
    # Analyze the reaction and determine the product
    # The reaction is a directed ortho-metalation followed by electrophilic quench.
    # The -CONEt2 group (position 1) directs to positions 2 and 6.
    # The -N(CH3)2 group (position 3) directs to positions 2 and 4.
    # The synergistic effect makes position 2 the most acidic and reactive site.
    # A methyl group is therefore added at position 2.
    final_product = "N,N-diethyl-2-methyl-3-dimethylaminobenzamide"
    
    # Print the reaction summary
    print("--- Reaction Summary ---")
    print(f"Starting Material: {starting_material}")
    print(f"Step 1 Reagents: {reagents_step1}")
    print(f"Step 2 Reagents: {reagents_step2}")
    print("-" * 25)
    
    # Print the explanation of how the product is formed
    print("\n--- Rationale ---")
    print("1. Directed Ortho-Metalation: The strong base, sec-BuLi, coordinated by TMEDA, removes the most acidic aromatic proton.")
    print("   - Both the amide group (at position 1) and the dimethylamino group (at position 3) direct the base.")
    print("   - The proton at position 2 is ortho to both groups, so it is selectively removed to form an aryllithium intermediate.")
    print("2. Electrophilic Quench: The aryllithium intermediate attacks the methyl group of methyl iodide.")
    print("   - A methyl group is added to position 2 of the benzene ring.")
    
    # Print the final product name
    print("\n--- Final Product ---")
    print("The compound obtained is:")
    print(final_product)

# Execute the function to display the result
chemical_reaction_synthesis()