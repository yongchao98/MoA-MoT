def solve_chemistry_problem():
    """
    Identifies and prints the reagents A and B for the given chemical reaction scheme.
    """

    reagent_A_name = "Hydrazine"
    reagent_A_formula = "H2N-NH2"

    reagent_B_name = "Propanal (Propionaldehyde)"
    reagent_B_formula = "CH3CH2CHO"

    print("The required reagents for the reaction scheme are:")
    print("-" * 45)
    
    # Description for the conversion of compound 1 to 2
    print("For the reaction from compound 1 to compound 2:")
    print(f"Reagent A is: {reagent_A_name} (Formula: {reagent_A_formula})")
    
    print("") # for spacing
    
    # Description for the conversion of compound 2 to 3
    print("For the reaction from compound 2 to compound 3:")
    print(f"Reagent B is: {reagent_B_name} (Formula: {reagent_B_formula})")
    print("-" * 45)

solve_chemistry_problem()