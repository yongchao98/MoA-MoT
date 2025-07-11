import sys

def solve_welding_parameters():
    """
    This script determines the appropriate TIG welding parameters for an Inconel 718
    turbine blade repair based on engineering principles.
    """

    # Given parameters from the problem description
    material = "Inconel 718"
    component_type = "Thin-section turbine blade tip"
    arc_gap_mm = 6.0

    # Answer choices
    options = {
        'A': {'current': 100.0, 'voltage': 7.5},
        'B': {'current': 17.5, 'voltage': 7.5},
        'C': {'current': 100.0, 'voltage': 15.0},
        'D': {'current': 150.0, 'voltage': 3.0},
        'E': {'current': 25.0, 'voltage': 17.5},
        'F': {'current': 80.0, 'voltage': 10.0}
    }

    print("### Analysis of Welding Parameters ###")
    print(f"Component: {component_type} made of {material}")
    print(f"Specified Arc Gap: {arc_gap_mm} mm\n")

    # --- Logic ---
    # Rule 1: For thin, heat-sensitive superalloys, a low current is required
    # to prevent overheating, distortion, and burn-through.
    # We define "low current" as being less than 50 A.
    low_current_threshold = 50.0
    print(f"Rule 1: Current should be low (e.g., < {low_current_threshold} A) for this delicate part.")

    # Rule 2: Voltage is proportional to arc length. A long arc requires a high voltage.
    # A 6mm arc is very long for TIG, so voltage should be high.
    # We define "high voltage" as being greater than 12 V for such a long arc.
    high_voltage_threshold = 12.0
    print(f"Rule 2: Voltage should be high (e.g., > {high_voltage_threshold} V) to sustain a long {arc_gap_mm} mm arc.\n")
    
    # Find the best option matching the criteria
    best_option_key = None
    for key, params in options.items():
        current = params['current']
        voltage = params['voltage']
        
        # Check if the option satisfies our engineering rules
        if current < low_current_threshold and voltage > high_voltage_threshold:
            best_option_key = key
            break

    if best_option_key:
        chosen_params = options[best_option_key]
        print("### Conclusion ###")
        print(f"Option '{best_option_key}' is the only choice with both low current and high voltage.")
        print("\nThe recommended parameters for the welding procedure are:")
        print(f"Current = {chosen_params['current']} A")
        print(f"Voltage = {chosen_params['voltage']} V")
        # Print the final answer in the required format
        print(f"\n<<<{best_option_key}>>>")
    else:
        print("Could not find a suitable option based on the defined engineering rules.")

solve_welding_parameters()