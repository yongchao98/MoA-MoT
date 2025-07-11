import pandas as pd

def solve_welding_parameters():
    """
    Analyzes welding parameters to find the most suitable option for TIG welding an Inconel 718 turbine blade.
    """
    # Given parameters
    travel_speed_mms = 0.5  # mm/s
    arc_gap_mm = 6.0        # mm

    # Answer Choices
    options = {
        'A': {'current_A': 100, 'voltage_V': 7.5},
        'B': {'current_A': 17.5, 'voltage_V': 7.5},
        'C': {'current_A': 100, 'voltage_V': 15},
        'D': {'current_A': 150, 'voltage_V': 3},
        'E': {'current_A': 25, 'voltage_V': 17.5},
        'F': {'current_A': 80, 'voltage_V': 10}
    }

    print("--- Step 1: Evaluating Voltage based on Arc Gap (6 mm) ---")
    print("A 6 mm arc gap is large for TIG welding and requires a relatively high voltage to sustain the arc.")
    print("Voltages below ~12V are generally considered too low for a stable 6 mm arc.\n")

    plausible_options = {}
    for option, params in options.items():
        voltage = params['voltage_V']
        # Rule of thumb: Voltages for a 6mm arc are typically > 12V.
        if voltage > 12:
            plausible_options[option] = params
            print(f"Option {option} ({voltage} V) has a plausible voltage.")
        else:
            print(f"Option {option} ({voltage} V) has a low/implausible voltage for a 6 mm arc gap.")

    print("\n--- Step 2: Calculating Heat Input for Plausible Options ---")
    print(f"Heat Input (J/mm) = (Voltage * Current) / Travel Speed ({travel_speed_mms} mm/s)")
    print("For Inconel 718, a low heat input (typically < 1000 J/mm) is crucial to prevent defects.\n")

    results = []
    for option, params in plausible_options.items():
        current = params['current_A']
        voltage = params['voltage_V']
        heat_input = (voltage * current) / travel_speed_mms
        results.append({
            'Option': option,
            'Current (A)': current,
            'Voltage (V)': voltage,
            'Heat Input (J/mm)': round(heat_input, 2)
        })

    results_df = pd.DataFrame(results)
    print(results_df.to_string(index=False))

    print("\n--- Step 3: Conclusion ---")
    print("Comparing the plausible options:")
    print("- Option C results in a very high heat input (3000.0 J/mm), which is unsuitable for Inconel 718 repair.")
    print("- Option E results in a controlled, low heat input (875.0 J/mm), which is ideal for this application.")
    print("\nTherefore, Option E provides the most appropriate parameters.")

    final_choice = 'E'
    final_params = options[final_choice]
    current = final_params['current_A']
    voltage = final_params['voltage_V']

    print("\nFinal Selected Parameters:")
    print(f"The recommended current is {current} A and the recommended voltage is {voltage} V.")

solve_welding_parameters()