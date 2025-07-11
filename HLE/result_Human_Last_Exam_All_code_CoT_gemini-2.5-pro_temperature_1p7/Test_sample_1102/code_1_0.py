import sys

def solve_welding_parameters():
    """
    This function determines the most appropriate welding parameters based on engineering principles
    and prints the reasoning and the result.
    """
    
    # Given parameters
    material = "Inconel 718"
    part = "Aeroengine compressor blade tip"
    process = "Manual TIG (GTAW)"
    arc_gap_mm = 6
    travel_speed_mm_s = 0.5

    # Best choice from the provided options
    chosen_current_A = 25
    chosen_voltage_V = 17.5
    
    print("Analysis for selecting the correct welding parameters:")
    print("1. Effect of Arc Gap on Voltage:")
    print(f"   - A long arc gap of {arc_gap_mm} mm requires a relatively high voltage to maintain a stable arc.")
    print("   - Based on TIG welding physics, the required voltage is estimated to be in the 15-20V range.")
    print(f"   - The voltage of {chosen_voltage_V} V is well-suited for this requirement.")
    print("\n2. Effect of Material and Speed on Current:")
    print(f"   - The part is a thin, heat-sensitive {material} blade tip.")
    print(f"   - The travel speed is very slow ({travel_speed_mm_s} mm/s), which increases heat input per unit length.")
    print("   - To prevent overheating and damage, a low welding current (amperage) is essential.")
    print(f"   - The current of {chosen_current_A} A is appropriately low for this delicate build-up operation.")
    
    print("\n" + "="*50)
    print("Conclusion: The recommended parameters for the welding procedure are:")
    # Final output as requested, printing each number in the final statement.
    print(f"Current = {chosen_current_A} A")
    print(f"Voltage = {chosen_voltage_V} V")
    print("="*50)

# Execute the function to solve the problem
solve_welding_parameters()