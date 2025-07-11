def solve_electric_field():
    """
    This function prints the derived formulas for the electric field in each region.
    """
    
    # The derived formula for the electric field in Region 1 (0 < phi < pi/2)
    E1_formula = "E_1 = (2 * sigma_2 * V_0) / (r * pi * (sigma_1 + sigma_2)) * i_phi"
    
    # The derived formula for the electric field in Region 2 (pi/2 < phi < pi)
    E2_formula = "E_2 = (2 * sigma_1 * V_0) / (r * pi * (sigma_1 + sigma_2)) * i_phi"

    print("The electric field in Region 1 is:")
    print(E1_formula)
    print("\nThe electric field in Region 2 is:")
    print(E2_formula)
    print("\nThese expressions match option C.")

solve_electric_field()