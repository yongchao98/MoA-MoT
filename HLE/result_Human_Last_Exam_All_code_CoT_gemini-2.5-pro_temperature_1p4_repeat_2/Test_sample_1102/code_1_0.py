import sys

def solve_welding_parameters():
    """
    Analyzes TIG welding parameters to find the most suitable setting for
    repairing an Inconel 718 turbine blade.
    """

    # --- Given Parameters ---
    travel_speed_mm_per_min = 30.0
    arc_gap_mm = 6.0
    material = "Inconel 718"
    component_type = "thin aerofoil tip"

    # --- Answer Choices: (Current in Amps, Voltage in Volts) ---
    choices = {
        'A': {'current': 100.0, 'voltage': 7.5},
        'B': {'current': 17.5, 'voltage': 7.5},
        'C': {'current': 100.0, 'voltage': 15.0},
        'D': {'current': 150.0, 'voltage': 3.0},
        'E': {'current': 25.0, 'voltage': 17.5},
        'F': {'current': 80.0, 'voltage': 10.0}
    }

    print("--- Step 1: Evaluating Voltage Plausibility ---")
    print(f"For a TIG arc gap of {arc_gap_mm} mm, a stable arc requires a sufficiently high voltage.")
    print("A typical voltage is > 15V. Let's assess each option:")
    
    plausible_choices = []
    for choice, params in choices.items():
        voltage = params['voltage']
        if voltage >= 15:
            status = "Plausible"
            plausible_choices.append(choice)
        else:
            status = "Unlikely (Too low to sustain a stable 6 mm arc)"
        print(f"Choice {choice}: Voltage = {voltage} V -> {status}")
    
    print("\n--- Step 2: Calculating Heat Input ---")
    print(f"Welding {material} on a {component_type} requires controlled (low) heat input to prevent defects.")
    print(f"Heat Input (J/mm) = (Voltage * Current * 60) / Travel Speed")
    
    best_choice = None
    min_plausible_heat_input = float('inf')

    for choice, params in choices.items():
        current = params['current']
        voltage = params['voltage']
        
        # Calculate heat input
        heat_input = (voltage * current * 60) / travel_speed_mm_per_min
        
        print(f"Choice {choice} ({current} A, {voltage} V): Heat Input = {heat_input:.1f} J/mm")

        # Determine the best choice among the plausible ones
        if choice in plausible_choices and heat_input < min_plausible_heat_input:
            min_plausible_heat_input = heat_input
            best_choice = choice
            
    print("\n--- Step 3: Final Conclusion ---")
    print(f"Choice C results in a very high heat input ({ (choices['C']['voltage'] * choices['C']['current'] * 60) / travel_speed_mm_per_min :.1f} J/mm), which is risky for this application.")
    print(f"Choice E provides a plausible voltage for the arc gap AND a much lower, more controllable heat input.")

    final_params = choices[best_choice]
    final_current = final_params['current']
    final_voltage = final_params['voltage']

    print(f"\nThe most suitable parameters are Current = {final_current} A and Voltage = {final_voltage} V.")
    print("The final calculation for the recommended heat input is:")
    print(f"({final_voltage} * {final_current} * 60) / {travel_speed_mm_per_min}")


solve_welding_parameters()
<<<E>>>