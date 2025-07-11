def calculate_displacement_field():
    """
    Calculates the displacement field in a dual-gate FET.

    The displacement field (D) is determined by the gate voltages and their
    respective capacitances per unit area. The formula used is:
    D = C_tg * V_tg + C_bg * V_bg
    where the transistor channel is grounded.
    """
    # --- User-definable parameters ---
    # Top gate capacitance per unit area (in Farads per square meter, F/m^2)
    Ctg = 1.5e-2

    # Back gate capacitance per unit area (in Farads per square meter, F/m^2)
    Cbg = 0.5e-2

    # Top gate voltage (in Volts, V)
    Vtg = 2.0

    # Back gate voltage (in Volts, V)
    Vbg = 5.0
    # -----------------------------------

    # Calculate the displacement field using the formula D = C_tg*V_tg + C_bg*V_bg
    displacement_field = Ctg * Vtg + Cbg * Vbg

    # Print the explanation and the step-by-step calculation
    print("The displacement field (D) is the total charge density induced in the channel.")
    print("The formula is: D = (Top Gate Capacitance) * (Top Gate Voltage) + (Back Gate Capacitance) * (Back Gate Voltage)")
    print("\nCalculating with the given values:")
    # Using f-string to format the output nicely, showing each number in the equation
    print(f"D = {Ctg:.2e} F/m^2 * {Vtg:.1f} V + {Cbg:.2e} F/m^2 * {Vbg:.1f} V")
    
    term1 = Ctg * Vtg
    term2 = Cbg * Vbg
    print(f"D = {term1:.2e} C/m^2 + {term2:.2e} C/m^2")

    # Print the final result. The unit of displacement field is Coulombs per square meter (C/m^2).
    print(f"\nThe total displacement field is: {displacement_field:.3e} C/m^2")

    return displacement_field

# Execute the function
final_D = calculate_displacement_field()
# The final answer is wrapped in <<<>>>
# print(f"<<<{final_D}>>>") # This would print the raw number, but let's make it formatted as in the output
print(f"\n<<<D = {final_D:.3e} C/m^2>>>")