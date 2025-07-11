import pandas as pd

def find_welding_parameters():
    """
    Analyzes potential TIG welding parameters for a turbine blade repair
    based on established engineering principles.
    """
    # Given parameters from the problem description
    material = "Inconel 718"
    process = "Manual TIG (GTAW) root pass build-up"
    arc_gap_mm = 6.0
    travel_speed_mm_s = 0.5

    # Answer choices
    options = [
        {'Choice': 'A', 'Current (A)': 100, 'Voltage (V)': 7.5},
        {'Choice': 'B', 'Current (A)': 17.5, 'Voltage (V)': 7.5},
        {'Choice': 'C', 'Current (A)': 100, 'Voltage (V)': 15},
        {'Choice': 'D', 'Current (A)': 150, 'Voltage (V)': 3},
        {'Choice': 'E', 'Current (A)': 25, 'Voltage (V)': 17.5},
        {'Choice': 'F', 'Current (A)': 80, 'Voltage (V)': 10}
    ]

    # --- Engineering Logic ---
    # Heuristic 1: Voltage vs. Arc Gap
    # A large arc gap of 6mm requires a high voltage to sustain the arc.
    # Typical TIG voltage is 10-20V. For 6mm, it should be > 14V.
    min_voltage_for_gap = 14.0

    # Heuristic 2: Current vs. Application
    # A root pass on a thin Inconel blade tip at a very slow travel speed
    # requires low current to avoid burn-through and overheating.
    max_current_for_application = 40.0

    print("Evaluating welding parameters based on the following criteria:")
    print(f"- Process: {process} on {material}")
    print(f"- Arc Gap: {arc_gap_mm} mm --> requires Voltage > {min_voltage_for_gap} V")
    print(f"- Application: Thin part, slow speed --> requires Current < {max_current_for_application} A\n")

    plausible_options = []
    for option in options:
        current = option['Current (A)']
        voltage = option['Voltage (V)']
        
        is_plausible = True
        reasons = []

        if voltage <= min_voltage_for_gap:
            is_plausible = False
            reasons.append(f"Voltage {voltage}V is too low for a {arc_gap_mm}mm arc gap.")
        
        if current >= max_current_for_application:
            is_plausible = False
            reasons.append(f"Current {current}A is too high for a delicate root pass.")

        if is_plausible:
            plausible_options.append(option)

    if len(plausible_options) == 1:
        best_option = plausible_options[0]
        print("Analysis complete. Found the most suitable parameters:")
        print(f"Choice {best_option['Choice']} is the only option that satisfies both criteria.")
        print("\n--- Recommended Welding Procedure Specification ---")
        # Final output of the "equation" as requested
        print(f"Current = {best_option['Current (A)']} A")
        print(f"Voltage = {best_option['Voltage (V)']} V")
        print("-------------------------------------------------")

    else:
        print("Could not determine a single best option based on the heuristics.")
        print("Plausible options found:", plausible_options)


if __name__ == "__main__":
    find_welding_parameters()