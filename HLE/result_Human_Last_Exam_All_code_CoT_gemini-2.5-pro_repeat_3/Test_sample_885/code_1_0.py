def solve_chemistry_problem():
    """
    This script provides the solution to the chemical retrosynthesis problem.
    It identifies the starting material in a Robinson annulation reaction.
    """

    starting_material_name = "ethyl 4-methyl-2-oxocyclohexanecarboxylate"
    
    # Chemical formulas for the balanced equation
    starting_material_formula = {"C": 10, "H": 16, "O": 2}
    mvk_formula = {"C": 4, "H": 6, "O": 1}
    product_formula = {"C": 14, "H": 20, "O": 3}
    water_formula = {"H": 2, "O": 1}

    print(f"The starting compound is: {starting_material_name}")
    print("\nThe overall reaction equation is:")
    
    # Printing the equation with all numbers as requested
    equation = (
        f"C{starting_material_formula['C']}H{starting_material_formula['H']}O{starting_material_formula['O']} + "
        f"C{mvk_formula['C']}H{mvk_formula['H']}O{mvk_formula['O']} -> "
        f"C{product_formula['C']}H{product_formula['H']}O{product_formula['O']} + "
        f"H{water_formula['H']}O{water_formula['O']}"
    )
    print(equation)

solve_chemistry_problem()