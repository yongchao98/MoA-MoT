import sys

def calculate_displacement_field(V_tg, V_bg, C_tg, C_bg):
    """
    Calculates the effective displacement field in a dual-gated FET.

    Args:
        V_tg (float): Top gate voltage in Volts (V).
        V_bg (float): Back gate voltage in Volts (V).
        C_tg (float): Top gate capacitance per unit area in Farads per square meter (F/m^2).
        C_bg (float): Back gate capacitance per unit area in Farads per square meter (F/m^2).

    Returns:
        float: The total displacement field in Coulombs per square meter (C/m^2).
    """
    # The displacement field is the sum of contributions from the top and bottom gates.
    # D = C_tg * V_tg + C_bg * V_bg
    D_total = C_tg * V_tg + C_bg * V_bg
    return D_total

def main():
    """
    Main function to execute the calculation and print the results.
    """
    # --- User-configurable example values ---
    # Top gate voltage in Volts
    V_tg = 1.0
    # Back gate voltage in Volts
    V_bg = 0.5
    # Top gate capacitance per area in F/m^2 (e.g., for ~2nm SiO2)
    C_tg = 1.72e-2
    # Bottom gate capacitance per area in F/m^2 (e.g., for ~25nm SiO2)
    C_bg = 0.138e-2
    # --- End of configuration ---

    print("Field-Effect Transistor Displacement Field Calculation\n")
    print("Physical Principle:")
    print("The effective displacement field (D) in the channel is the sum of the contributions from the top and back gates.")
    print("The formula is: D = C_tg * V_tg + C_bg * V_bg\n")

    print("Given values:")
    print(f"  Top Gate Voltage (V_tg)    = {V_tg} V")
    print(f"  Back Gate Voltage (V_bg)   = {V_bg} V")
    print(f"  Top Gate Capacitance (C_tg)  = {C_tg:.4g} F/m^2")
    print(f"  Back Gate Capacitance (C_bg) = {C_bg:.4g} F/m^2\n")

    # Calculate the displacement field
    D = calculate_displacement_field(V_tg, V_bg, C_tg, C_bg)

    # Print the final equation with all numbers
    print("Calculation:")
    # The final print statement fulfills the requirement to output each number in the final equation.
    print(f"D = {C_tg} * {V_tg} + {C_bg} * {V_bg}")

    # Print the result
    # Units of D are (F/m^2) * V = (C/V/m^2) * V = C/m^2
    print(f"\nThe resulting displacement field (D) is: {D:.4g} C/m^2")

    # Final answer in the required format for automated checking
    sys.stdout.write(f"\n<<<{D}>>>")

if __name__ == "__main__":
    main()