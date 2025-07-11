def solve_displacement_field():
    """
    Calculates and prints the symbolic equation for the displacement field in a dual-gate FET.
    """
    # Define the symbolic variables. We use strings to represent them for the final printout.
    Ctg = "C_tg"  # Top gate capacitance per area
    Cbg = "C_bg"  # Back gate capacitance per area
    Vtg = "V_tg"  # Top gate voltage
    Vbg = "V_bg"  # Back gate voltage

    # The displacement field (D) is the superposition of the fields created by the top and back gates.
    # The field from the back gate is D_bg = C_bg * V_bg.
    # The field from the top gate is D_tg = -C_tg * V_tg (assuming positive direction is from back to top).
    # The total displacement field is D = D_bg + D_tg.

    # Print the final equation, showing each term clearly as requested.
    print("The displacement field (D) through the transistor is given by the superposition of the fields from the top and back gates.")
    print("The final equation is:")
    print(f"D = {Cbg} * {Vbg} - {Ctg} * {Vtg}")

solve_displacement_field()
