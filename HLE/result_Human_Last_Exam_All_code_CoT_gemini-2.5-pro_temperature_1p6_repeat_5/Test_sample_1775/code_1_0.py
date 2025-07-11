import sys

def calculate_displacement_field(Ctg, Vtg, Cbg, Vbg):
    """
    Calculates the displacement field in a dual-gate FET.

    Args:
        Ctg (float): Top gate capacitance per unit area (in F/m^2).
        Vtg (float): Top gate voltage (in V).
        Cbg (float): Back gate capacitance per unit area (in F/m^2).
        Vbg (float): Back gate voltage (in V).

    Returns:
        float: The calculated displacement field (in C/m^2).
    """
    # The displacement field D is the sum of the contributions from the top and back gates.
    # D = Ctg * Vtg + Cbg * Vbg
    displacement_field = Ctg * Vtg + Cbg * Vbg
    return displacement_field

def main():
    """
    Main function to execute the calculation and print the results.
    """
    # --- Example Values ---
    # Top gate capacitance per area in Farads per square meter (F/m^2)
    # Example: 1.5 uF/cm^2 = 1.5e-6 F / (1e-2 m)^2 = 1.5e-2 F/m^2
    Ctg = 1.5e-2
    # Top gate voltage in Volts (V)
    Vtg = 1.0
    # Back gate capacitance per area in Farads per square meter (F/m^2)
    # Example: 0.5 uF/cm^2 = 0.5e-6 F / (1e-2 m)^2 = 0.5e-2 F/m^2
    Cbg = 0.5e-2
    # Back gate voltage in Volts (V)
    Vbg = -2.0

    # Calculate the displacement field
    D = calculate_displacement_field(Ctg, Vtg, Cbg, Vbg)

    # Print the explanation and the final equation with values
    print("The displacement field (D) through the transistor is the sum of the fields induced by the top and back gates.")
    print("The formula is: D = Ctg * Vtg + Cbg * Vbg")
    print("\nUsing the provided values:")
    # sys.stdout.write is used to prevent adding an extra space before the plus sign
    sys.stdout.write(f"D = ({Ctg:.2e} F/m^2) * ({Vtg:.2f} V) + ({Cbg:.2e} F/m^2) * ({Vbg:.2f} V)\n")
    print("\nResult:")
    # The unit of displacement field is Coulombs per square meter (C/m^2)
    print(f"D = {D:.2e} C/m^2")
    
    # This line is for the final answer extraction and is not part of the standard output.
    # It prints the raw numerical result for parsing.
    print(f"Final_Answer_For_Extraction: {D}")


if __name__ == '__main__':
    main()
