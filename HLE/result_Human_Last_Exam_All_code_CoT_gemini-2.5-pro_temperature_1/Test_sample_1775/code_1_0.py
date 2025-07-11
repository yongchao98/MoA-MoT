def calculate_displacement_field(Ctg, Vtg, Cbg, Vbg):
    """
    Calculates the net displacement field in a dual-gate FET channel.

    Args:
        Ctg (float): Top gate capacitance per unit area (in F/m^2).
        Vtg (float): Top gate voltage (in V).
        Cbg (float): Back gate capacitance per unit area (in F/m^2).
        Vbg (float): Back gate voltage (in V).
    """
    # The displacement field is the sum of the fields induced by each gate.
    # D = D_top + D_back = (Ctg * Vtg) + (Cbg * Vbg)
    D = (Ctg * Vtg) + (Cbg * Vbg)

    # Print the final equation with the numbers plugged in
    print("The formula for the displacement field (D) is: D = (Ctg * Vtg) + (Cbg * Vbg)")
    print("Plugging in the values:")
    # Using the 'f-string' feature to format the output string
    # The format specifier 'e' is used for scientific notation
    print(f"D = ({Ctg:.2e} F/m^2 * {Vtg:.2f} V) + ({Cbg:.2e} F/m^2 * {Vbg:.2f} V)")
    
    # Print the intermediate and final results
    term1 = Ctg * Vtg
    term2 = Cbg * Vbg
    print(f"D = {term1:.2e} C/m^2 + {term2:.2e} C/m^2")
    print(f"D = {D:.2e} C/m^2")
    
    return D

# --- Example Usage ---
# You can change these values to match your specific problem.
# Typical capacitance values are in the range of uF/cm^2, which is 10^-2 F/m^2.
# Let's use some example values.
example_Ctg = 1.5e-2  # F/m^2, equivalent to 1.5 uF/cm^2
example_Vtg = 1.0     # Volts
example_Cbg = 0.5e-2  # F/m^2, equivalent to 0.5 uF/cm^2
example_Vbg = -2.0    # Volts

# Calculate and print the result
final_D = calculate_displacement_field(example_Ctg, example_Vtg, example_Cbg, example_Vbg)

# To provide a single numerical answer for the example above
# final_D = (1.5e-2 * 1.0) + (0.5e-2 * -2.0) = 1.5e-2 - 1.0e-2 = 0.5e-2
# <<<0.005>>>