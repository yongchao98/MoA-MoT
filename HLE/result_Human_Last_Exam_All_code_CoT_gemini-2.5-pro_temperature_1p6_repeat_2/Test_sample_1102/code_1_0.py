def calculate_welding_parameters():
    """
    Calculates and explains the optimal welding parameters for the given scenario.
    """
    # Given parameters
    material = "Inconel 718"
    process = "Manual TIG (GTAW) build-up"
    travel_speed_mms = 0.5  # mm/s
    arc_gap_mm = 6  # mm

    # --- Parameter Analysis ---
    # Voltage is primarily determined by arc gap. A long 6mm gap requires a high voltage.
    # A typical estimate is V_arc = V_drop + (E * L).
    # Using V_drop=10V and E=1.5V/mm: 10V + (1.5V/mm * 6mm) = 19V.
    # Voltages of 15V or 17.5V are plausible.

    # Current determines heat input. For a slow build-up on a sensitive superalloy,
    # a low current is needed to control the weld pool size and prevent overheating.
    # Heat Input (J/mm) = (Voltage * Current) / Travel Speed

    # Chosen parameters from Option E
    chosen_voltage = 17.5  # Volts
    chosen_current = 25.0    # Amps

    # --- Calculation ---
    # We select Option E (25 A, 17.5 V) as it provides a high enough voltage for the arc gap
    # and a low enough current for a controlled build-up at a slow travel speed.
    heat_input = (chosen_voltage * chosen_current) / travel_speed_mms

    print("Welding Parameter Justification:")
    print(f"Process: {process} on {material}")
    print(f"Arc Gap: {arc_gap_mm} mm, Travel Speed: {travel_speed_mms} mm/s")
    print("\nA long arc gap requires a high voltage to be sustained.")
    print("A slow travel speed on a sensitive superalloy requires a low current to control heat input.")
    print("\nBased on this, the chosen parameters are:")
    print(f"Current: {chosen_current} A")
    print(f"Voltage: {chosen_voltage} V")

    print("\nTo confirm this, we calculate the linear heat input:")
    print("Formula: Heat Input = (Voltage * Current) / Travel Speed")
    # Outputting the numbers in the final equation
    print(f"Calculation: ({chosen_voltage} V * {chosen_current} A) / {travel_speed_mms} mm/s = {heat_input:.2f} J/mm")
    print("\nThis heat input is controlled and suitable for a precise build-up, making these parameters the correct choice.")

calculate_welding_parameters()