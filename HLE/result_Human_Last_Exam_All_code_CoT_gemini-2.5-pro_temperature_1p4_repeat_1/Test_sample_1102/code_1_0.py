import sys

def solve_welding_parameters():
    """
    Analyzes welding parameters for a specific TIG repair process and outputs the recommended settings.
    """
    # Known parameters from the problem description
    material = "Inconel 718"
    application = "Aeroengine turbine blade tip repair"
    travel_speed_mms = 0.5  # Travel speed in mm/s
    arc_gap_mm = 6.0        # Arc gap in mm

    # --- Reasoning ---
    # 1. Voltage vs. Arc Gap: A 6 mm arc gap is large for TIG welding and requires a
    #    relatively high voltage to be sustained. Typical TIG voltages are above 10V.
    #    This makes choices with voltages of 3V, 7.5V, or 10V highly unlikely.
    #    We are left with plausible options C (15V) and E (17.5V).

    # 2. Current vs. Application: The part is a thin turbine blade made of Inconel 718,
    #    requiring precise heat control. The travel speed (0.5 mm/s) is very slow, which
    #    increases heat concentration.
    #    - A high current like 100A (Option C) would introduce excessive heat, likely
    #      damaging or melting through the delicate blade tip.
    #    - A low current like 25A (Option E) is much more appropriate for such a
    #      delicate repair, allowing for controlled material build-up without
    #      overheating the base metal.

    # 3. Conclusion: The combination of low current (25A) and a necessarily high
    #    voltage (17.5V) to support the long arc is the most logical and safe choice.
    
    # Selected parameters based on the reasoning
    chosen_current_A = 25.0
    chosen_voltage_V = 17.5

    # Calculate the Heat Input in Joules per millimeter
    # Heat Input (J/mm) = (Voltage * Current) / Travel Speed (mm/s)
    heat_input_J_per_mm = (chosen_voltage_V * chosen_current_A) / travel_speed_mms

    # --- Final Output ---
    print("Based on the analysis of the welding requirements, the recommended parameters are:")
    print("-" * 40)
    print(f"Selected Current: {int(chosen_current_A)} A")
    print(f"Selected Voltage: {chosen_voltage_V} V")
    print("-" * 40)
    
    print("\nThis choice provides a low current suitable for a delicate Inconel blade,")
    print("paired with a voltage high enough to sustain the specified 6 mm arc gap.")
    
    print("\nThe heat input for this procedure is calculated as follows:")
    # Outputting each number in the final equation as requested
    print(f"Heat Input Equation: ({chosen_voltage_V} V * {int(chosen_current_A)} A) / {travel_speed_mms} mm/s")
    print(f"Result: {heat_input_J_per_mm:.1f} J/mm")

# Execute the function
solve_welding_parameters()