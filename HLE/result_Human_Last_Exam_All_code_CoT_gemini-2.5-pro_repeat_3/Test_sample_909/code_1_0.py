def solve_electric_field():
    """
    This function prints the derived expressions for the electric field
    in the two regions of the resistor.
    """
    # The derived expression for the electric field in Region 1 (0 < phi < pi/2)
    E1_expression = "E_1 = (2 * sigma_2 * V_0) / (r * pi * (sigma_1 + sigma_2)) * i_phi"

    # The derived expression for the electric field in Region 2 (pi/2 < phi < pi)
    E2_expression = "E_2 = (2 * sigma_1 * V_0) / (r * pi * (sigma_1 + sigma_2)) * i_phi"

    print("The electric field in the first region (0 < phi < pi/2) is:")
    print(E1_expression)
    print("\nThe electric field in the second region (pi/2 < phi < pi) is:")
    print(E2_expression)
    print("\nThis corresponds to Answer Choice C.")

solve_electric_field()