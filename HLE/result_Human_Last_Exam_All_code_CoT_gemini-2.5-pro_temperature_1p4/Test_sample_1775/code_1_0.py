def calculate_displacement_field(C_tg, V_tg, C_bg, V_bg):
    """
    Calculates the displacement field in a dual-gate transistor.

    Args:
        C_tg (float): Top gate capacitance per unit area (e.g., in F/cm^2).
        V_tg (float): Top gate voltage (in V).
        C_bg (float): Bottom gate capacitance per unit area (e.g., in F/cm^2).
        V_bg (float): Bottom gate voltage (in V).

    Returns:
        float: The calculated displacement field (e.g., in C/cm^2).
    """
    return (C_tg * V_tg) + (C_bg * V_bg)

def main():
    """
    Main function to demonstrate the calculation with example values.
    """
    # Example values for demonstration
    # Capacitance per area in Farads per square centimeter (F/cm^2)
    C_tg = 1.725e-6  # Represents ~2 nm SiO2 gate dielectric
    C_bg = 0.8625e-6 # Represents ~4 nm SiO2 gate dielectric

    # Gate voltages in Volts (V)
    V_tg = 1.5
    V_bg = -2.0

    # Calculate the displacement field
    displacement_field = calculate_displacement_field(C_tg, V_tg, C_bg, V_bg)
    
    # Components of the calculation
    D_tg = C_tg * V_tg
    D_bg = C_bg * V_bg

    print("The formula for the displacement field (D) is:")
    print("D = (C_tg * V_tg) + (C_bg * V_bg)")
    print("\nPlugging in the example values:")
    # Print the full equation with numbers
    print(f"D = ({C_tg} F/cm^2 * {V_tg} V) + ({C_bg} F/cm^2 * {V_bg} V)")
    print(f"D = {D_tg:.4g} C/cm^2 + {D_bg:.4g} C/cm^2")
    print(f"\nFinal Result:")
    # Use scientific notation for the final result for clarity
    print(f"D = {displacement_field:.4g} C/cm^2")

if __name__ == "__main__":
    main()