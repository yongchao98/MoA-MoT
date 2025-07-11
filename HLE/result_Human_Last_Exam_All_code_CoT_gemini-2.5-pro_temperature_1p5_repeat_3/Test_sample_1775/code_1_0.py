def solve_displacement_field():
    """
    This function prints the equation for the displacement field in a dual-gate FET.
    """
    # Symbolic variable names
    C_tg = "C_tg"  # Top gate capacitance per area
    V_tg = "V_tg"  # Top gate voltage
    C_bg = "C_bg"  # Back gate capacitance per area
    V_bg = "V_bg"  # Back gate voltage

    # The displacement field D in the transistor is determined by the total charge
    # density induced in the channel by both gates. This is the sum of the
    # charge densities induced by each gate individually.
    
    # Equation: D = (Charge from Top Gate) + (Charge from Back Gate)
    # D = (C_tg * V_tg) + (C_bg * V_bg)
    
    print("The equation for the displacement field (D) is:")
    print(f"D = {C_tg} * {V_tg} + {C_bg} * {V_bg}")

solve_displacement_field()