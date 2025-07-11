def solve_reaction():
    """
    Determines the IUPAC name of the product from the reaction of
    1,3-dibromo-2-iodobenzene and excess phenyl magnesium bromide.
    """

    # 1. Define the reactants and conditions
    reactant_aromatic = "1,3-dibromo-2-iodobenzene"
    reactant_grignard = "excess phenyl magnesium bromide (Ph-MgBr)"
    conditions = "reflux in THF, then aqueous work-up"

    print("Step 1: Identify the reactants and relative reactivity.")
    print(f"The reaction is between {reactant_aromatic} and {reactant_grignard}.")
    print("A Grignard reagent (Ph-MgBr) reacts with an aryl halide.")
    print("The reactivity of the C-Halogen bonds is: C-I > C-Br > C-Cl.")
    print("-" * 20)

    # 2. Predict the first substitution
    print("Step 2: First substitution reaction.")
    print("The most reactive halogen, iodine (I) at position 2, is replaced by a phenyl group.")
    intermediate_product_1 = "1,3-dibromo-2-phenylbenzene"
    print(f"Product after first substitution: {intermediate_product_1}")
    print("-" * 20)

    # 3. Predict the subsequent substitutions
    print("Step 3: Reaction with excess Grignard reagent.")
    print("Because excess phenyl magnesium bromide is used under reflux (heat), the less reactive bromine (Br) atoms are also replaced.")
    final_product_structure = "1,2,3-triphenylbenzene"
    print("The bromine atoms at positions 1 and 3 are replaced by phenyl groups.")
    print("-" * 20)

    # 4. Determine the final IUPAC name
    final_iupac_name = "1,2,3-triphenylbenzene"
    print("Final Step: Determine the IUPAC name of the final product.")
    print("The structure has a central benzene ring with three phenyl groups at positions 1, 2, and 3.")
    print("\nThe IUPAC name of the final product is:")
    print(final_iupac_name)

solve_reaction()