import sys

def solve_welding_parameters():
    """
    Calculates and analyzes welding parameters to determine the best choice
    for TIG welding a Rolls-Royce Trent 1000 compressor blade tip.
    """

    # --- Given Parameters ---
    travel_speed_mm_per_min = 30.0  # Given as 0.5 mm/s, which is 0.5 * 60 = 30 mm/min
    arc_gap_mm = 6.0
    material = "Inconel 718"
    application = "Thin compressor blade tip root pass"

    # --- Analysis ---
    print(f"Analyzing welding parameters for {material} on a {application}.")
    print(f"Given Travel Speed: {travel_speed_mm_per_min} mm/min")
    print(f"Given Arc Gap: {arc_gap_mm} mm\n")

    # Step 1: Estimate the required arc voltage based on the arc gap.
    # A rule of thumb for TIG welding is V_arc = 10 + (0.5 * Arc_Length_mm).
    # This provides a baseline for required voltage.
    voltage_constant = 10.0
    voltage_factor = 0.5
    estimated_voltage = voltage_constant + (voltage_factor * arc_gap_mm)
    
    print("--- Step 1: Voltage Analysis ---")
    print("A long arc gap requires a high voltage to maintain a stable arc.")
    print("Using the rule-of-thumb V = 10 + (0.5 * Arc Gap):")
    print(f"Estimated Voltage = {voltage_constant} + ({voltage_factor} * {arc_gap_mm}) = {estimated_voltage:.1f} V")
    print("This suggests that choices with high voltages (like 15 V or 17.5 V) are most plausible.\n")

    # Step 2: Analyze heat input for critical candidates.
    # Inconel 718 requires low heat input (typically < 1.5 kJ/mm) to prevent defects.
    # Heat Input (J/mm) = (Voltage * Current * 60) / Travel_Speed (mm/min)
    print("--- Step 2: Heat Input Analysis ---")
    print(f"Welding {material} requires low heat input to prevent defects like cracking.")
    
    # Candidate C: 100 A and 15 V
    current_C = 100.0
    voltage_C = 15.0
    heat_input_C = (voltage_C * current_C * 60) / travel_speed_mm_per_min

    print("\nAnalysis of Candidate C (100 A, 15 V):")
    print(f"Heat Input = ({voltage_C} V * {current_C} A * 60) / {travel_speed_mm_per_min} mm/min = {heat_input_C / 1000:.3f} kJ/mm")
    print("Result: This heat input is too high for a delicate Inconel 718 repair and risks damaging the part.")

    # Candidate E: 25 A and 17.5 V
    current_E = 25.0
    voltage_E = 17.5
    heat_input_E = (voltage_E * current_E * 60) / travel_speed_mm_per_min
    
    print("\nAnalysis of Candidate E (25 A, 17.5 V):")
    print(f"Heat Input = ({voltage_E} V * {current_E} A * 60) / {travel_speed_mm_per_min} mm/min = {heat_input_E / 1000:.3f} kJ/mm")
    print("Result: This heat input is low and controlled, which is ideal for Inconel 718.")
    
    # --- Conclusion ---
    print("\n--- Final Conclusion ---")
    print("Option E (25 A and 17.5 V) is the only choice that satisfies both critical requirements:")
    print(f"1. The voltage ({voltage_E} V) is appropriately high for the long {arc_gap_mm} mm arc gap.")
    print(f"2. The resulting heat input ({heat_input_E / 1000:.3f} kJ/mm) is low and suitable for the heat-sensitive {material}.")

solve_welding_parameters()