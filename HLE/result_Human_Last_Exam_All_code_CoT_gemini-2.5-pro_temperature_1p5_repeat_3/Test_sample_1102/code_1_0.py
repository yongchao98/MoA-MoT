def calculate_welding_parameters():
    """
    Analyzes and calculates a key welding parameter to determine the correct settings
    for TIG welding a thin Inconel 718 turbine blade tip.

    The analysis is based on the following logic:
    1. Material: Inconel 718 blade tip is thin and heat-sensitive, requiring LOW current.
    2. Arc Gap: The specified 6 mm arc gap is very large for TIG, requiring HIGH voltage.
    3. Conclusion: The combination of low current and high voltage is needed. From the choices,
       25 A and 17.5 V is the most suitable combination.

    This script calculates the resulting heat input for these parameters.
    """
    # Given parameters
    travel_speed_mm_per_min = 30.0  # 0.5 mm/s * 60 s/min

    # Chosen parameters from the most suitable answer choice (E)
    current_amps = 25.0
    voltage_volts = 17.5

    # Welding efficiency for TIG (GTAW) is typically around 0.6 to 0.8
    # For this calculation, we use a common value of 0.7
    efficiency = 0.7

    # Calculate Heat Input (HI) in Joules per millimeter (J/mm)
    # HI = (Voltage * Current * Efficiency * 60) / Travel Speed (mm/min)
    heat_input = (voltage_volts * current_amps * efficiency * 60) / travel_speed_mm_per_min

    print("Based on welding principles for thin Inconel 718 and a large arc gap:")
    print(f"Suggested Current (I): {current_amps} A")
    print(f"Suggested Voltage (V): {voltage_volts} V")
    print("-" * 30)
    print("This provides a controlled heat input suitable for the repair.")
    print("Calculating the Heat Input:")
    # The prompt requires outputting each number in the final equation.
    print(f"Equation: (V * I * Efficiency * 60) / S_mm_per_min")
    print(f"Calculation: ({voltage_volts} V * {current_amps} A * {efficiency} * 60) / {travel_speed_mm_per_min} mm/min")
    print(f"Resulting Heat Input: {heat_input:.2f} J/mm")

calculate_welding_parameters()
