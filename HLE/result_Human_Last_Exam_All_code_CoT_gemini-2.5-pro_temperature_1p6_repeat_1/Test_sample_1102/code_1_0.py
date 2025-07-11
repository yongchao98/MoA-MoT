def solve_welding_parameters():
    """
    Analyzes TIG welding parameters for a turbine blade repair and identifies the most suitable option.
    """
    # Given parameters from the problem description
    travel_speed_mmpmin = 30.0  # mm/min (equivalent to 0.5 mm/s)
    arc_gap_mm = 6.0

    # Parameters from the most suitable answer choice (E)
    current_amps = 25.0
    voltage_volts = 17.5

    # Explanation
    print("Based on the analysis of the welding requirements for a thin Inconel 718 turbine blade:")
    print(f"- A large arc gap of {arc_gap_mm} mm requires a relatively high voltage to sustain the arc.")
    print(f"- A slow travel speed of {travel_speed_mmpmin} mm/min requires a low current to prevent excessive heat input and damage.")
    print("\nOption E (25 A and 17.5 V) is the best fit as it satisfies both conditions.")
    print("-" * 60)

    # Calculate and display the heat input
    heat_input_j_per_mm = (voltage_volts * current_amps * 60) / travel_speed_mmpmin

    print("Calculating the Heat Input for the chosen parameters:")
    print(f"Heat Input = ({voltage_volts} V * {current_amps} A * 60) / {travel_speed_mmpmin} mm/min")
    print(f"Resulting Heat Input = {heat_input_j_per_mm:.0f} J/mm")
    print("\nThis heat input value is appropriate for a controlled build-up on a delicate superalloy part.")

    print("\nTherefore, the recommended parameters for the welding procedure are:")
    print(f"Current: {int(current_amps)} A")
    print(f"Voltage: {voltage_volts} V")


solve_welding_parameters()