def identify_reaction_product():
    """
    This function traces a multi-step organic synthesis reaction to identify the final product.
    It prints the name of the chemical at each step of the reaction.
    """
    # Step 1: Friedel-Crafts Acylation
    # Benzene + Propanoyl chloride --(AlCl3)--> Intermediate-1
    intermediate_1 = "1-phenylpropan-1-one"
    print("Step 1: Friedel-Crafts Acylation")
    print(f"The first intermediate is: {intermediate_1}\n")

    # Step 2: Electrophilic Aromatic Bromination
    # Intermediate-1 --(Br2/FeBr3)--> Intermediate-2
    # The acyl group is a meta-director.
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    print("Step 2: Electrophilic Aromatic Bromination")
    print(f"The second intermediate is: {intermediate_2}\n")

    # Step 3: Catalytic Hydrogenation
    # Intermediate-2 --(H2/Pd)--> Intermediate-3
    # The ketone is reduced to a methylene group.
    intermediate_3 = "1-bromo-3-propylbenzene"
    print("Step 3: Catalytic Hydrogenation (Ketone Reduction)")
    print(f"The third intermediate is: {intermediate_3}\n")

    # Step 4: Free-Radical Benzylic Bromination
    # Intermediate-3 --(NBS, (PhCO2)2)--> Final Product
    # Bromination occurs at the benzylic position.
    final_product = "1-bromo-1-(3-bromophenyl)propane"
    print("Step 4: Free-Radical Benzylic Bromination")
    print(f"The final product is: {final_product}")

identify_reaction_product()