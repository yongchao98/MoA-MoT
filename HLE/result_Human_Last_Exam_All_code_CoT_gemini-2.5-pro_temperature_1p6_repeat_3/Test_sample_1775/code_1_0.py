def solve_displacement_field():
    """
    This function prints the derived equation for the displacement field in the FET.
    It uses placeholder variable names to construct and display the final formula.
    """
    # Define the variables as strings for the equation
    Ctg = "Ctg"  # Top gate capacitance per area
    Vtg = "Vtg"  # Top gate voltage
    Cbg = "Cbg"  # Back gate capacitance per area
    Vbg = "Vbg"  # Back gate voltage

    # The equation for the displacement field D
    # D = Ctg*Vtg + Cbg*Vbg
    
    print("The displacement field (D) is given by the equation:")
    # We print each component of the equation as requested.
    print(f"D = ({Ctg} * {Vtg}) + ({Cbg} * {Vbg})")

solve_displacement_field()