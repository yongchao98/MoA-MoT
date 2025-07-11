import sys

def solve_welding_parameters():
    """
    Analyzes TIG welding parameters to find the suitable choice for repairing
    an Inconel 718 turbine blade tip.
    """
    # --- Given Parameters ---
    material = "Inconel 718"
    application = "Tip repair of a thin aeroengine blade"
    arc_gap_mm = 6.0
    travel_speed_mm_per_min = 30.0

    # --- Answer Choices ---
    # Format: (Choice_Letter, Current_Amps, Voltage_Volts)
    choices = [
        ('A', 100, 7.5),
        ('B', 17.5, 7.5),
        ('C', 100, 15),
        ('D', 150, 3),
        ('E', 25, 17.5),
        ('F', 80, 10)
    ]

    print("Analyzing welding parameter choices based on physical principles...")
    print(f"Material: {material}, Application: {application}")
    print(f"Arc Gap: {arc_gap_mm} mm, Travel Speed: {travel_speed_mm_per_min} mm/min\n")

    # --- Evaluation Logic ---
    # Rule 1: Voltage for a 6mm arc must be high.
    # A TIG arc voltage is approx. V_drop + V_gradient * length.
    # (e.g., 10V + 1.5V/mm * 6mm = 19V). So, we expect V > 15V.
    min_voltage = 15.0
    max_voltage = 22.0

    # Rule 2: Current for a delicate build-up on a thin Inconel part
    # at slow speed should be low to control heat input.
    min_current = 15.0
    max_current = 50.0

    best_choice = None
    best_choice_letter = ''

    for choice, current, voltage in choices:
        print(f"--- Evaluating Choice {choice}: {current} A, {voltage} V ---")
        
        # Check Voltage
        if voltage >= min_voltage and voltage <= max_voltage:
            voltage_ok = True
            print(f"  [PASS] Voltage ({voltage} V) is plausible for a long ({arc_gap_mm} mm) arc.")
        else:
            voltage_ok = False
            print(f"  [FAIL] Voltage ({voltage} V) is too low to sustain a {arc_gap_mm} mm arc. Expected >{min_voltage} V.")

        # Check Current
        if current >= min_current and current <= max_current:
            current_ok = True
            print(f"  [PASS] Current ({current} A) is appropriate for delicate repair work.")
        else:
            current_ok = False
            print(f"  [FAIL] Current ({current} A) is outside the expected range for precise tip build-up.")

        if voltage_ok and current_ok:
            best_choice = (current, voltage)
            best_choice_letter = choice
        print("")


    if best_choice:
        final_current, final_voltage = best_choice
        print(f"Conclusion: Choice {best_choice_letter} is the only option that satisfies both voltage and current requirements.")
        print(f"The most suitable parameters are {final_current} A and {final_voltage} V.\n")
        
        # As requested, showing the numbers in a final equation (Heat Input)
        heat_input = (final_voltage * final_current * 60) / travel_speed_mm_per_min
        print("This corresponds to a controlled heat input, calculated as follows:")
        
        # Outputting each number in the final equation as requested
        print(f"Heat Input = ({final_voltage} V * {final_current} A * 60) / {travel_speed_mm_per_min} mm/min = {heat_input:.2f} J/mm")
        
        # Hide the final answer from the regular output to put it at the end
        sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf8', closefd=False) # reopen stdout
        sys.stdout.write(f"\n<<<{best_choice_letter}>>>")
    else:
        print("No suitable choice found based on the defined criteria.")

solve_welding_parameters()