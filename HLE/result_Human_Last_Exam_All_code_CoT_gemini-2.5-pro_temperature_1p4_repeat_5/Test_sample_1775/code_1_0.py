def solve_displacement_field():
    """
    This function prints the symbolic formula for the displacement field
    in a dual-gate field-effect transistor.
    """
    # Symbolic representation of the variables
    C_bg = "C_bg"  # Back gate capacitance per area
    V_bg = "V_bg"  # Back gate voltage
    C_tg = "C_tg"  # Top gate capacitance per area
    V_tg = "V_tg"  # Top gate voltage

    # The displacement field D is the superposition of the fields from the back and top gates.
    # Field from back gate: C_bg * V_bg
    # Field from top gate: -C_tg * V_tg (negative sign due to direction)
    # Total D = (C_bg * V_bg) - (C_tg * V_tg)

    print("The displacement field (D) is given by the equation:")
    # Print the final equation showing all the symbolic components
    print(f"D = ({C_bg} * {V_bg}) - ({C_tg} * {V_tg})")

solve_displacement_field()