def solve_reaction():
    """
    This function determines the product of the reaction between 1,3-dibromo-2-iodobenzene
    and excess phenyl magnesium bromide and prints its IUPAC name.
    """
    # Reactants and Conditions
    starting_material = "1,3-dibromo-2-iodobenzene"
    reagent = "Phenyl magnesium bromide (PhMgBr)"
    conditions = "Excess reagent and reflux"

    # Analysis of reactivity
    # In reactions with Grignard reagents, the C-I bond is much more reactive than the C-Br bond.
    # The reaction proceeds by replacing the most reactive halogen first.
    # Given the excess reagent and heat (reflux), the reaction is driven to completion,
    # replacing all halogen atoms with phenyl groups.

    # Step-wise substitution
    # Step 1: Iodine is replaced by a phenyl group.
    intermediate_1 = "1,3-dibromo-2-phenylbenzene"
    
    # Step 2: One of the bromine atoms is replaced by a phenyl group.
    intermediate_2 = "1-bromo-2,3-diphenylbenzene"

    # Step 3: The final bromine atom is replaced by a phenyl group.
    final_product = "1,2,3-triphenylbenzene"

    print(f"The reaction starts with {starting_material}.")
    print("Under excess phenyl magnesium bromide and reflux conditions, all halogens are substituted.")
    print("The order of substitution is based on halogen reactivity (I > Br).")
    print(f"1. Substitution of Iodine leads to: {intermediate_1}")
    print(f"2. Substitution of a Bromine leads to: {intermediate_2}")
    print(f"3. Substitution of the final Bromine leads to the final product.")
    
    print("\nThe IUPAC name of the final product is:")
    
    # Printing the final name, including the numbers as requested.
    name_parts = final_product.split('-')
    numbers = name_parts[0]
    substituent = name_parts[1]
    parent_chain = name_parts[2]

    print(f"{numbers}-{substituent}-{parent_chain}")

solve_reaction()