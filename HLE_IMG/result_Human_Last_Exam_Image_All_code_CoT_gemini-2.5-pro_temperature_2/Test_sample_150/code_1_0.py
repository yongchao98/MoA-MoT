def solve_reaction_pathway():
    """
    Identifies the products at each step of the given chemical reaction sequence.
    """
    start_molecule_1 = "Benzene"
    start_molecule_2 = "Propanoyl chloride"
    
    intermediate_1 = "Propiophenone (1-phenylpropan-1-one)"
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    intermediate_3 = "1-(3-bromophenyl)propane"
    final_product = "1-bromo-1-(3-bromophenyl)propane"

    print(f"Step 1: Friedel-Crafts acylation of {start_molecule_1} with {start_molecule_2}.")
    print(f"Intermediate-1 is: {intermediate_1}\n")
    
    print(f"Step 2: Electrophilic bromination of Intermediate-1.")
    print(f"Intermediate-2 is: {intermediate_2}\n")

    print(f"Step 3: Hydrogenation (reduction) of the ketone in Intermediate-2.")
    print(f"Intermediate-3 is: {intermediate_3}\n")
    
    print(f"Step 4: Benzylic bromination of Intermediate-3 with NBS.")
    print(f"The Final Product is: {final_product}\n")

solve_reaction_pathway()