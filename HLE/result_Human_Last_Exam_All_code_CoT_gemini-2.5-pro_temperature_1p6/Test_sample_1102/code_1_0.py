import sys

def calculate_welding_parameters():
    """
    This function analyzes welding parameters for a specific repair scenario
    and determines the most appropriate current and voltage settings.
    """
    # Given parameters
    material = "Inconel 718"
    process = "Manual TIG (GTAW) root pass"
    component = "Rolls-Royce Trent 1000 compressor blade tip"
    arc_gap_mm = 6.0
    travel_speed_mm_per_min = 30.0

    print("Analyzing Welding Procedure for Turbine Blade Repair\n")
    print(f"Component: {component} ({material})")
    print(f"Process: {process}")
    print(f"Specified Arc Gap: {arc_gap_mm} mm (long arc)")
    print(f"Specified Travel Speed: {travel_speed_mm_per_min} mm/min (slow)\n")

    print("Step 1: Evaluate Voltage Requirement")
    print(f"A long TIG arc of {arc_gap_mm} mm requires a relatively high voltage to be sustained.")
    print("Options with low voltages (e.g., < 12 V) are less likely.")
    print("Options C (15 V) and E (17.5 V) are the most plausible candidates based on voltage.\n")

    print("Step 2: Calculate and Compare Heat Input for Plausible Options")
    print("Heat Input is critical for preventing defects in thin-section Inconel 718.")
    print("Formula: Heat Input (J/mm) = (Voltage * Current * 60) / Travel Speed (mm/min)\n")

    # Option C: 100 A and 15 V
    current_C = 100.0
    voltage_C = 15.0
    heat_input_C = (voltage_C * current_C * 60) / travel_speed_mm_per_min

    print(f"Analysis of Option C ({current_C} A, {voltage_C} V):")
    print(f"Heat Input = ({voltage_C} V * {current_C} A * 60) / {travel_speed_mm_per_min} mm/min = {heat_input_C:.0f} J/mm")
    print("Result: This heat input is very high and would likely cause overheating, distortion, or cracking in a delicate blade tip.\n")

    # Option E: 25 A and 17.5 V
    current_E = 25.0
    voltage_E = 17.5
    heat_input_E = (voltage_E * current_E * 60) / travel_speed_mm_per_min

    print(f"Analysis of Option E ({current_E} A, {voltage_E} V):")
    print(f"Heat Input = ({voltage_E} V * {current_E} A * 60) / {travel_speed_mm_per_min} mm/min = {heat_input_E:.0f} J/mm")
    print("Result: This is a much lower, controlled heat input appropriate for sensitive aerospace repairs.\n")

    print("Step 3: Conclusion")
    print("The combination of high voltage and low current in Option E is a coherent strategy for this specific repair:")
    print(f"- The high voltage ({voltage_E} V) is necessary to maintain the specified long ({arc_gap_mm} mm) arc.")
    print(f"- The low current ({current_E} A) provides a low heat input to prevent damage to the thin blade material.")
    print("This low-power, wide-arc technique allows for stable material build-up without excessive penetration.\n")
    
    print("Final Answer: The recommended parameters are from Option E.")
    print("The final equation for the recommended heat input is:")
    print(f"Heat Input = {voltage_E} V * {current_E} A * 60 / {travel_speed_mm_per_min} mm/min = {heat_input_E:.0f} J/mm")

if __name__ == '__main__':
    calculate_welding_parameters()