def calculate_displacement_field(Ctg, Vtg, Cbg, Vbg):
    """
    Calculates the net displacement field in a dual-gate transistor.

    Args:
        Ctg (float): Top gate capacitance per unit area (in F/m^2).
        Vtg (float): Top gate voltage (in V).
        Cbg (float): Back gate capacitance per unit area (in F/m^2).
        Vbg (float): Back gate voltage (in V).
    """

    # The formula for the net displacement field (D) is derived from the superposition
    # of the fields from the top and back gates.
    # We define the positive direction as pointing from the back gate to the top gate.
    # D = (Displacement from Back Gate) - (Displacement from Top Gate)
    # D = C_bg * V_bg - C_tg * V_tg

    D_net = Cbg * Vbg - Ctg * Vtg

    # Print the explanation and the final equation with the numbers plugged in.
    print("The formula for the net displacement field (D) is: D = C_bg * V_bg - C_tg * V_tg")
    print("\nPlugging in the given values:")
    # The f-string below constructs the equation string with the actual numbers.
    # The signs are handled naturally. A negative number will appear as '- C*V'.
    # To make it clear like "A * B - C * D", we handle the sign of the second term.
    # This is for display purposes only.
    if Vtg >= 0:
        print(f"D = ({Cbg}) * ({Vbg}) - ({Ctg}) * ({Vtg})")
    else:
        # If Vtg is negative, the minus sign is already there.
        print(f"D = ({Cbg}) * ({Vbg}) - ({Ctg}) * ({Vtg})")
        
    print(f"\nCalculated Result:")
    print(f"D = {Cbg * Vbg} - {Ctg * Vtg}")
    print(f"D = {D_net} (in C/m^2)")


if __name__ == '__main__':
    # Since no values were provided in the problem, we will use example values
    # to demonstrate the calculation as requested.
    # All units are SI units.
    example_Ctg = 1.7e-2  # F/m^2 (equivalent to a ~2 nm thick SiO2 layer)
    example_Vtg = 1.0     # Volts
    example_Cbg = 0.5e-2  # F/m^2 (equivalent to a ~7 nm thick SiO2 layer)
    example_Vbg = -0.5    # Volts

    calculate_displacement_field(example_Ctg, example_Vtg, example_Cbg, example_Vbg)