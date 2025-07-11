def solve_reaction_product_name():
    """
    This script determines the IUPAC name of the product from the reaction of
    1,3-dibromo-2-iodobenzene with excess phenyl magnesium bromide.
    """

    # Step 1: Define the reactants and conditions.
    # Substrate: 1,3-dibromo-2-iodobenzene
    # Reagent: Excess phenyl magnesium bromide (PhMgBr)
    # Conditions: Reflux in THF, followed by aqueous work-up.

    # Step 2: Determine the reaction pathway.
    # Given the excess Grignard reagent and heat (reflux), a complete
    # substitution of all halogen atoms by phenyl groups is the most likely outcome.
    # The reaction proceeds via sequential substitution, with reactivity C-I > C-Br.

    # Reaction sequence:
    # 1. 1,3-dibromo-2-iodobenzene + PhMgBr -> 1,3-dibromo-2-phenylbenzene + MgBrI
    # 2. 1,3-dibromo-2-phenylbenzene + PhMgBr -> 3-bromo-1,2-diphenylbenzene + MgBr2
    # 3. 3-bromo-1,2-diphenylbenzene + PhMgBr -> 1,2,3-triphenylbenzene + MgBr2

    # Step 3: Identify the final product and its IUPAC name.
    # The final product has a benzene ring with three phenyl groups attached
    # at the 1, 2, and 3 positions.
    final_iupac_name = "1,2,3-Triphenylbenzene"
    
    print("The IUPAC name of the major product is:")
    print(final_iupac_name)
    
    # As requested, output each number from the final name's locants.
    print("\nThe locant numbers in the name are:")
    locants = ['1', '2', '3']
    for number in locants:
        print(number)

# Execute the function to get the answer.
solve_reaction_product_name()