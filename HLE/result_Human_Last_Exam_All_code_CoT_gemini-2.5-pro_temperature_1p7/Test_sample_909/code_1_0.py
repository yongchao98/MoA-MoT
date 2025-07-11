def solve_electric_field():
    """
    This function prints the derived expressions for the electric field
    in the two regions of the cylindrical resistor.
    The derivation is based on solving Laplace's equation with the given
    boundary conditions.
    """
    
    # The derived expressions for the electric fields in the two regions are:
    E1_expression = "(2 * sigma_2 * V_0) / (r * pi * (sigma_1 + sigma_2))"
    E2_expression = "(2 * sigma_1 * V_0) / (r * pi * (sigma_1 + sigma_2))"
    
    print("The problem is solved by applying Laplace's equation and boundary conditions.")
    print("The resulting electric field in each region is purely in the phi-direction (i_phi).")
    
    # --- Output for Region 1 ---
    print("\n--- Electric Field in Region 1 (0 < phi < pi/2) ---")
    print(f"E_1 = ({E1_expression}) * i_phi")
    
    print("\nBreaking down the expression for E_1:")
    print("  Numerator: 2 * sigma_2 * V_0")
    print("  Denominator: r * pi * (sigma_1 + sigma_2)")
    print("  The terms are:")
    print("    '2' and 'pi': Mathematical constants")
    print("    'V_0': Applied DC voltage")
    print("    'r': Radial distance from the center")
    print("    'sigma_1', 'sigma_2': Ohmic conductivities of the two regions")

    # --- Output for Region 2 ---
    print("\n--- Electric Field in Region 2 (pi/2 < phi < pi) ---")
    print(f"E_2 = ({E2_expression}) * i_phi")

    print("\nBreaking down the expression for E_2:")
    print("  Numerator: 2 * sigma_1 * V_0")
    print("  Denominator: r * pi * (sigma_1 + sigma_2)")
    print("  The terms are the same as listed above.")
    
    print("\nThese results correspond to Answer Choice C.")

solve_electric_field()