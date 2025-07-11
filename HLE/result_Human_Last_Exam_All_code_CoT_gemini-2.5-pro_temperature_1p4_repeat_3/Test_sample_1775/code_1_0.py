def solve_displacement_field():
    """
    This function calculates and prints the formula for the displacement field
    in a dual-gate FET.
    """
    # Symbolic representation of variables
    D = "D"
    C_tg = "C_tg"
    V_tg = "V_tg"
    C_bg = "C_bg"
    V_bg = "V_bg"

    print("The net displacement field (D) through the transistor is the superposition of the fields from the top and back gates.")
    print("We define the positive direction as pointing from the back gate to the top gate.")
    print("\nThe field from the back gate points in the positive direction.")
    print("The field from the top gate points in the negative direction.")
    print("\nThe final equation for the displacement field is:")
    # Using f-string to format and print the equation with variables
    print(f"{D} = ({C_bg} * {V_bg}) - ({C_tg} * {V_tg})")

    print("\nWhere:")
    print(f"  {D} = Net displacement field")
    print(f"  {C_bg} = Back gate capacitance per area")
    print(f"  {V_bg} = Back gate voltage")
    print(f"  {C_tg} = Top gate capacitance per area")
    print(f"  {V_tg} = Top gate voltage")

# Execute the function
solve_displacement_field()