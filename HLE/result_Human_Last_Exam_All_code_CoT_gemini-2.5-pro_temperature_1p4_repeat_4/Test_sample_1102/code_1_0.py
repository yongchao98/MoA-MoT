def find_welding_parameters():
    """
    Analyzes welding parameters for TIG welding of an Inconel 718 turbine blade
    and selects the most appropriate settings from a list of choices.
    """
    # Given parameters from the problem description
    material = "Inconel 718"
    process = "Manual TIG (GTAW)"
    travel_speed_mm_per_sec = 0.5
    travel_speed_mm_per_min = travel_speed_mm_per_sec * 60
    arc_gap_mm = 6

    print("### Welding Parameter Analysis ###\n")
    print(f"Material: {material}")
    print(f"Process: {process}")
    print(f"Application: Turbine blade tip repair (root pass build-up)")
    print(f"Travel Speed: {travel_speed_mm_per_sec} mm/s ({travel_speed_mm_per_min} mm/min)")
    print(f"Arc Gap: {arc_gap_mm} mm\n")

    print("### Evaluation of Choices ###")
    print("1. Voltage Requirement: A 6 mm arc gap is large for TIG welding and requires a relatively high voltage to maintain a stable arc. Voltages below ~12V are highly unlikely to sustain such an arc. This eliminates choices with very low voltages (3V, 7.5V).")
    print("2. Current Requirement: For a root pass build-up on an Inconel turbine blade tip, the current must be sufficient for fusion without causing excessive heat input, which could lead to defects. A very low current (<50A) may result in a 'cold' weld, while a very high current (>120A) risks burn-through.")
    print("\nBased on this, let's assess the options:")
    print(" - A (100A, 7.5V), B (17.5A, 7.5V), D (150A, 3V): Eliminated due to unrealistically low voltage for a 6mm arc.")
    print(" - E (25A, 17.5V): Voltage is plausible, but 25A is likely too low for proper fusion and build-up on this component.")
    print(" - F (80A, 10V): While 80A is a reasonable current, 10V is too low for a stable 6mm arc gap.")
    print(" - C (100A, 15V): This option provides a realistic voltage (15V) to sustain a 6mm arc and a sufficient current (100A) to ensure good fusion for the build-up without being excessive.\n")

    print("### Conclusion ###")
    print("The most suitable parameters are 100 Amps and 15 Volts.")

    # Chosen parameters from option C
    chosen_current_A = 100  # Amps
    chosen_voltage_V = 15    # Volts

    # Calculate the heat input using the chosen parameters
    # Heat Input (kJ/mm) = (Voltage * Amperage * 0.06) / Travel Speed (mm/min)
    heat_input = (chosen_voltage_V * chosen_current_A * 0.06) / travel_speed_mm_per_min

    print("\nThe recommended welding procedure specification should list:")
    print(f"Current: {chosen_current_A} A")
    print(f"Voltage: {chosen_voltage_V} V\n")

    print("This corresponds to a calculated heat input for the process.")
    # The final print statement shows each number in the equation as requested
    print(f"Heat Input Calculation: ({chosen_voltage_V} V * {chosen_current_A} A * 0.06) / {travel_speed_mm_per_min} mm/min = {heat_input:.2f} kJ/mm")

# Execute the function to get the answer
find_welding_parameters()