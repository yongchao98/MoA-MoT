def solve_electric_field():
    """
    This function prints the derived symbolic expressions for the electric field
    in the two regions of the cylindrical resistor.
    """
    # The electric field expressions are derived from first principles
    # (Laplace's equation with boundary conditions).
    
    # E_vector = E_magnitude * phi_hat
    # where phi_hat is the unit vector in the azimuthal direction.
    
    # Expression for the electric field in Region 1 (0 < phi < pi/2)
    E1_expression = "E1 = (2 * sigma_2 * V_0) / (r * pi * (sigma_1 + sigma_2)) * i_phi"
    
    # Expression for the electric field in Region 2 (pi/2 < phi < pi)
    E2_expression = "E2 = (2 * sigma_1 * V_0) / (r * pi * (sigma_1 + sigma_2)) * i_phi"
    
    print("The electric field in Region 1 is:")
    print(E1_expression)
    print("\nThe electric field in Region 2 is:")
    print(E2_expression)

solve_electric_field()