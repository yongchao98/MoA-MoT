import sys

def find_welding_parameters():
    """
    This script determines the appropriate TIG welding parameters for repairing an
    Inconel 718 turbine blade tip based on engineering principles.
    """

    # --- Given Information ---
    material = "Inconel 718"
    application = "Tip repair of a thin aeroengine blade"
    arc_gap_mm = 6.0
    travel_speed_mms = 0.5

    print("Step 1: Analyze the required Arc Voltage based on the Arc Gap.")
    print(f"The specified arc gap is {arc_gap_mm} mm, which is exceptionally long for TIG welding.")
    print("Arc voltage is directly proportional to the arc gap length.")
    print("A common approximation is: V_arc = V_drop + (E_field * Length)")
    print("For TIG, this is roughly V ~= 10-12V + (1-2 V/mm * Length).")
    estimated_voltage_min = 10 + (1 * arc_gap_mm)
    estimated_voltage_max = 12 + (2 * arc_gap_mm)
    print(f"This implies a required voltage in the range of {estimated_voltage_min:.1f}V to {estimated_voltage_max:.1f}V to sustain the arc.")
    print("Based on this, we can eliminate choices with very low voltages (< 15V).\n")

    print("Step 2: Analyze the required Current based on the application.")
    print(f"The application is a delicate repair on a thin '{material}' blade tip.")
    print("This requires low heat input to prevent warping, burn-through, and material damage.")
    print(f"The travel speed is very slow ({travel_speed_mms} mm/s), which already increases heat input per unit length.")
    print("Therefore, a low welding current is necessary for a controlled, stable repair.\n")

    print("Step 3: Evaluate the most plausible option.")
    print("Comparing the likely choices:")
    print(" - Option C (100 A, 15 V): The voltage is plausible, but 100 A is an extremely high current for this delicate task and would likely cause damage.")
    print(" - Option E (25 A, 17.5 V): The voltage of 17.5 V is ideal for the long arc gap. The current of 25 A is appropriately low for precise control on a thin section.")
    print("\nConclusion: The combination of a low current for control and a high voltage to sustain the long arc is the correct choice.\n")

    # --- Final Recommended Parameters ---
    chosen_current_A = 25
    chosen_voltage_V = 17.5

    print("Final Recommended Parameters for the Welding Procedure Specification:")
    print(f"Current = {chosen_current_A} A")
    print(f"Voltage = {chosen_voltage_V} V")

find_welding_parameters()
# Appending the final answer choice as requested by the format.
# The double slash is just to make it a comment in python, the final output will be the answer marker.
sys.stdout.write("<<<E>>>\n")