import sys

def solve_welding_parameters():
    """
    Analyzes TIG welding parameters to select the correct specification.
    """
    # Given parameters from the problem description
    arc_gap_mm = 6
    material = "Inconel 718"
    application = "Build-up repair root pass"

    # Answer choices
    choices = {
        "A": {"current_A": 100, "voltage_V": 7.5},
        "B": {"current_A": 17.5, "voltage_V": 7.5},
        "C": {"current_A": 100, "voltage_V": 15},
        "D": {"current_A": 150, "voltage_V": 3},
        "E": {"current_A": 25, "voltage_V": 17.5},
        "F": {"current_A": 80, "voltage_V": 10},
    }

    print("### Step-by-Step Analysis ###")
    
    # --- Step 1: Analyze Voltage Requirement ---
    print(f"\n1. Analyzing Voltage based on Arc Gap ({arc_gap_mm} mm):")
    # Rule of thumb: V_arc ≈ 10 + 0.5 * L_mm
    estimated_voltage = 10 + (0.5 * arc_gap_mm)
    voltage_lower_bound = 12.0 # A reasonable minimum for a stable long arc
    print(f"   - A long arc gap of {arc_gap_mm} mm requires a relatively high voltage to remain stable.")
    print(f"   - A common estimation is: Voltage ≈ 10 + (0.5 * Arc Length) = 10 + (0.5 * {arc_gap_mm}) = {estimated_voltage} V.")
    print(f"   - Therefore, we will filter for options with a voltage of at least {voltage_lower_bound} V.")

    # --- Step 2: Analyze Current Requirement ---
    print(f"\n2. Analyzing Current for a '{application}':")
    current_lower_bound = 50 # Amps
    current_upper_bound = 120 # Amps
    print(f"   - The application requires sufficient current to melt the {material} and filler wire.")
    print(f"   - A very low current (< {current_lower_bound} A) would be insufficient for material build-up.")
    print(f"   - A very high current (> {current_upper_bound} A) risks burn-through on a blade tip.")
    print(f"   - A current between {current_lower_bound} A and {current_upper_bound} A is most suitable.")

    # --- Step 3: Evaluate and Select the Best Choice ---
    print("\n3. Evaluating Choices:")
    best_choice = None
    for key, params in choices.items():
        current = params['current_A']
        voltage = params['voltage_V']
        
        # Check if parameters are plausible
        is_voltage_ok = voltage >= voltage_lower_bound
        is_current_ok = current_lower_bound <= current <= current_upper_bound
        
        if is_voltage_ok and is_current_ok:
            best_choice = key
            break # Assume only one choice will fit the criteria

    # --- Step 4: Final Answer Output ---
    if best_choice:
        final_params = choices[best_choice]
        final_current = final_params['current_A']
        final_voltage = final_params['voltage_V']

        print("\n### Conclusion ###")
        print(f"The most suitable parameters are in Option {best_choice}.")
        print(f"The voltage of {final_voltage} V is appropriate for the {arc_gap_mm} mm arc gap.")
        print(f"The current of {final_current} A is suitable for the build-up repair application.")
        
        print("\nFinal selected numbers for the procedure:")
        print(f"Current = {final_current}")
        print(f"Voltage = {final_voltage}")
        
        # Output for the platform
        print(f"\n<<<{best_choice}>>>", file=sys.stderr)
    else:
        print("\nCould not determine a single best choice based on the criteria.")


solve_welding_parameters()