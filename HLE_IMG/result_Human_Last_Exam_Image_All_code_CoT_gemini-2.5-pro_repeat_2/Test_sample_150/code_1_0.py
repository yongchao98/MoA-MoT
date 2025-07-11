def identify_product():
    """
    This function tracks the transformations in the provided reaction scheme
    and identifies the final product.
    """
    # Step 1: Friedel-Crafts Acylation
    reactant_1 = "Benzene"
    reactant_2 = "Propanoyl chloride"
    intermediate_1 = "Propiophenone"
    print(f"Step 1: {reactant_1} reacts with {reactant_2} in a Friedel-Crafts acylation.")
    print(f"Intermediate-1 is: {intermediate_1}\n")

    # Step 2: Electrophilic Aromatic Bromination
    reagent_2 = "Br2/FeBr3"
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    print(f"Step 2: {intermediate_1} is brominated with {reagent_2}.")
    print("The acyl group is a meta-director.")
    print(f"Intermediate-2 is: {intermediate_2}\n")

    # Step 3: Catalytic Hydrogenation (Reduction)
    reagent_3 = "H2/Pd"
    intermediate_3 = "1-bromo-3-propylbenzene"
    print(f"Step 3: {intermediate_2} is reduced with {reagent_3}.")
    print("The ketone group is reduced to an alkyl group.")
    print(f"Intermediate-3 is: {intermediate_3}\n")

    # Step 4: Benzylic Bromination
    reagent_4 = "NBS"
    final_product = "1-bromo-1-(3-bromophenyl)propane"
    print(f"Step 4: {intermediate_3} undergoes benzylic bromination with {reagent_4}.")
    print("A bromine atom is added to the carbon directly attached to the benzene ring.")
    print(f"The final product is: {final_product}")

identify_product()