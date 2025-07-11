def get_iupac_name():
    """
    Determines the IUPAC name of the product from the reaction of
    1,3-dibromo-2-iodobenzene with excess phenyl magnesium bromide.
    """
    reactant = "1,3-dibromo-2-iodobenzene"
    reagent = "excess phenyl magnesium bromide"

    print(f"Reaction: {reactant} with {reagent}")
    print("------------------------------------------")

    # Step 1: Substitution of the most reactive halogen (Iodine)
    print("Step 1: The carbon-iodine bond is the weakest, so the iodine is the best leaving group.")
    print("The first phenyl group from the Grignard reagent substitutes the iodine at position 2.")
    intermediate_1 = "1,3-dibromo-2-phenylbenzene"
    print(f"Intermediate product: {intermediate_1}\n")

    # Step 2: Substitution of the first bromine
    print("Step 2: Since excess reagent is used under reflux (heating), the less reactive bromine atoms will also be substituted.")
    print("A second phenyl group substitutes one of the bromine atoms (e.g., at position 1).")
    intermediate_2 = "3-bromo-1,2-diphenylbenzene"
    print(f"Intermediate product: {intermediate_2}\n")

    # Step 3: Substitution of the final bromine
    print("Step 3: A third phenyl group substitutes the final remaining bromine atom at position 3.")
    final_product_description = "A benzene ring with three phenyl substituents at positions 1, 2, and 3."
    print(f"Final product: {final_product_description}\n")
    
    # Step 4: Aqueous work-up
    print("Step 4: The final aqueous work-up quenches any remaining Grignard reagent and does not change the organic product.\n")

    # Construct the IUPAC name
    print("Constructing the final IUPAC name:")
    locants = ["1", "2", "3"]
    prefix = "Triphenyl"
    parent = "benzene"
    
    # Print each number in the final name as requested
    print(f"The locants (positions) for the phenyl groups are: {', '.join(locants)}.")
    
    final_name = f"{locants[0]},{locants[1]},{locants[2]}-{prefix}{parent}"
    
    print("\n--- Final IUPAC Name ---")
    print(final_name)


get_iupac_name()
<<<1,2,3-Triphenylbenzene>>>