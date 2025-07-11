def calculate_displacement_field():
    """
    Calculates the displacement field in a dual-gate FET.

    The displacement field (D) through the transistor is the sum of the
    contributions from the top and back gates, equivalent to the total
    charge density induced in the grounded channel.
    The formula is: D = C_tg * V_tg + C_bg * V_bg
    """

    # --- Define Example Values ---
    # C_tg: Top gate capacitance per unit area (in Farads per square meter, F/m^2)
    C_tg = 2.4e-3
    # V_tg: Top gate voltage (in Volts, V)
    V_tg = 1.5
    # C_bg: Back gate capacitance per unit area (in F/m^2)
    C_bg = 6.0e-4
    # V_bg: Back gate voltage (in V)
    V_bg = 5.0

    # --- Calculation ---
    displacement_field = C_tg * V_tg + C_bg * V_bg
    term1 = C_tg * V_tg
    term2 = C_bg * V_bg

    # --- Output ---
    print("The displacement field (D) is calculated using the formula: D = C_tg * V_tg + C_bg * V_bg")
    print("\nGiven the example values:")
    print(f"C_tg = {C_tg} F/m^2")
    print(f"V_tg = {V_tg} V")
    print(f"C_bg = {C_bg} F/m^2")
    print(f"V_bg = {V_bg} V")
    print("\nThe calculation is:")
    # Print the equation with all the numbers, as requested
    print(f"D = ({C_tg}) * ({V_tg}) + ({C_bg}) * ({V_bg})")
    print(f"D = {term1} + {term2}")
    print(f"D = {displacement_field} C/m^2")

# Run the calculation and print the result
calculate_displacement_field()