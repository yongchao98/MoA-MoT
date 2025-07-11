def solve_displacement_field_equation():
    """
    This function prints the final symbolic equation for the displacement
    field (D) in a dual-gate FET.
    """
    # Symbolic variables for the final equation
    C_tg = "C_tg"
    V_tg = "V_tg"
    C_bg = "C_bg"
    V_bg = "V_bg"

    # The equation for the displacement field D is the sum of the
    # charge densities induced by the top and back gates.
    # D = (charge from top gate) + (charge from back gate)
    # D = (C_tg * V_tg) + (C_bg * V_bg)

    print("The displacement field (D) through the transistor is given by the following equation:")
    # We print each component of the final equation as requested.
    print(f"D = ({C_tg}) * ({V_tg}) + ({C_bg}) * ({V_bg})")

# Execute the function to display the result.
solve_displacement_field_equation()