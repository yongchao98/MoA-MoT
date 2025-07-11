def solve_welding_parameters():
    """
    Analyzes TIG welding parameters to select the best option for repairing an
    Inconel 718 turbine blade tip.
    """
    # --- Given Parameters ---
    travel_speed_mm_per_min = 30.0
    arc_gap_mm = 6.0
    material = "Inconel 718"
    application = "Turbine blade tip repair (build-up)"

    # --- Answer Choices ---
    choices = {
        'A': {'current_A': 100, 'voltage_V': 7.5},
        'B': {'current_A': 17.5, 'voltage_V': 7.5},
        'C': {'current_A': 100, 'voltage_V': 15},
        'D': {'current_A': 150, 'voltage_V': 3},
        'E': {'current_A': 25, 'voltage_V': 17.5},
        'F': {'current_A': 80, 'voltage_V': 10}
    }

    # --- Analysis ---
    print("### Welding Parameter Analysis ###")
    print(f"Objective: Select optimal TIG welding parameters for {application} on {material}.")
    print(f"Given Conditions: Travel Speed = {travel_speed_mm_per_min} mm/min, Arc Gap = {arc_gap_mm} mm.\n")

    print("--- Step 1: Voltage Evaluation ---")
    print(f"A TIG arc gap of {arc_gap_mm} mm is unusually long. A stable arc requires a relatively high voltage.")
    print("Options with very low voltage (< 12 V) are generally not feasible for this arc length.")
    print("This makes options A, B, D, and F technically unlikely.\n")
    
    print("--- Step 2: Heat Input Comparison for Plausible Options (C and E) ---")
    print(f"For a heat-sensitive material like {material} on a thin blade, low heat input is crucial to prevent defects.")
    print("Heat Input Formula: HI (J/mm) = (Voltage * Current * 60) / Travel Speed\n")

    # Analysis of Choice C
    c_params = choices['C']
    c_current, c_voltage = c_params['current_A'], c_params['voltage_V']
    c_heat_input = (c_voltage * c_current * 60) / travel_speed_mm_per_min
    print(f"Analysis of Option C ({c_current} A, {c_voltage} V):")
    print(f"Heat Input = ({c_voltage} V * {c_current} A * 60) / {travel_speed_mm_per_min} mm/min = {c_heat_input:.0f} J/mm")
    print("This represents a high heat input, posing a risk of overheating and defects for this application.\n")
    
    # Analysis of Choice E
    e_params = choices['E']
    e_current, e_voltage = e_params['current_A'], e_params['voltage_V']
    e_heat_input = (e_voltage * e_current * 60) / travel_speed_mm_per_min
    print(f"Analysis of Option E ({e_current} A, {e_voltage} V):")
    print(f"Heat Input = ({e_voltage} V * {e_current} A * 60) / {travel_speed_mm_per_min} mm/min = {e_heat_input:.0f} J/mm")
    print("This is a much lower, more controlled heat input suitable for precision build-up work.\n")

    print("--- Conclusion ---")
    print("Option E provides the best combination:")
    print(f"1. The voltage ({e_voltage} V) is appropriate for the long {arc_gap_mm} mm arc gap.")
    print(f"2. The low current and slow speed result in low heat input, ideal for a successful repair on Inconel 718.")
    
    print("\n### Final Recommended Parameters ###")
    print(f"Current: {e_current} A")
    print(f"Voltage: {e_voltage} V")

solve_welding_parameters()