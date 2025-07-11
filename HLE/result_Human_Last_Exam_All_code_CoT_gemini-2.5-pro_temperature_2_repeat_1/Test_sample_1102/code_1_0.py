import sys

def solve_welding_parameters():
    """
    Analyzes the provided TIG welding scenario to determine the appropriate current and voltage,
    then prints the reasoning and the selected answer.
    """

    # Given parameters
    material = "Inconel 718"
    application = "Aeroengine turbine blade tip repair (root pass)"
    travel_speed_mms = 0.5
    arc_gap_mm = 6

    # --- Analysis ---

    # 1. Current analysis based on heat input requirements
    # The application involves a thin part and slow travel speed, requiring low heat input.
    # High currents (>= 80A) are unsuitable. Low current is required.
    # Of the choices, 17.5A and 25A are low current options.
    chosen_current = 25  # Amps

    # 2. Voltage analysis based on arc gap
    # TIG voltage is primarily a function of arc length.
    # A simplified equation is V = V_fall + E_gradient * Arc_Length
    # For Argon, V_fall is ~10-12V and E_gradient is ~1-1.5 V/mm.
    # Estimated Voltage = 10V + (1.5 V/mm * 6 mm) = 19 V.
    # We need a voltage in this vicinity. 7.5V and 10V are too low for a 6mm arc. 17.5V is the most plausible value.
    chosen_voltage = 17.5  # Volts

    # --- Output ---
    print("Based on the principles of TIG (GTAW) welding for the specified repair:")
    print(f"\n1. Current (Amps): The part is a thin Inconel 718 blade tip. To prevent overheating, distortion, or burn-through, especially at a slow travel speed of {travel_speed_mms} mm/s, a low heat input is critical. A high current would be destructive. Therefore, a low current is required.")

    print(f"\n2. Voltage (Volts): The arc voltage is primarily dependent on the arc gap, which is specified as {arc_gap_mm} mm. A stable arc over this relatively large distance requires a sufficiently high voltage. A voltage around 17-20 V is expected. Lower voltages would likely fail to sustain the arc.")

    print("\nConclusion: The combination of a low current and a moderately high voltage is the only physically sound option for a successful repair under these conditions.")

    print("\n--- Final Equation ---")
    print("The selected parameters for the welding procedure specification are:")
    print(f"Current = {chosen_current} A")
    print(f"Voltage = {chosen_voltage} V")

solve_welding_parameters()
<<<E>>>