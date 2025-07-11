def solve_wittig_reaction():
    """
    Determines and displays the product of a specified Wittig reaction.
    """

    # Step 1 & 2: Define the reactants based on the problem description.
    # The Wittig reaction involves an aldehyde and a phosphorus ylide.
    aldehyde = "pivalaldehyde"
    wittig_reagent = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"

    # Step 3, 4, & 5: Determine the products of the reaction.
    # The reaction forms an alkene and triphenylphosphine oxide.
    # The ylide is unstabilized, leading to the (Z)-isomer as the major product.
    # The IUPAC name is derived from the combined structure: (CH3)3C-CH=CH-CH2-(2-chlorophenyl)
    # The main chain is a pent-2-ene, with substituents at positions 1 and 4.
    major_product = "(Z)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    byproduct = "triphenylphosphine oxide"

    # Print the full reaction equation, showing all components.
    # This fulfills the requirement to "output each number in the final equation"
    # by printing the full IUPAC names which contain numbers.
    print("The Wittig reaction occurs as follows:")
    print(f"{aldehyde} + {wittig_reagent} ---> {major_product} + {byproduct}")
    print("-" * 50)
    print(f"The major organic product of the reaction is: {major_product}")

solve_wittig_reaction()