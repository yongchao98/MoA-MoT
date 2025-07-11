import math

def calculate_reactive_power_compensation():
    """
    Calculates the total reactive power compensation required for an HVDC system.

    The calculation is based on the following assumptions:
    1. Reactive power compensation for harmonics is a percentage of the active power,
       derived from the given harmonic distortion levels.
    2. Reactive power compensation for voltage stability after a fault is also
       a percentage of the active power, derived from the given voltage drop percentage.
    """

    # Given values from the problem description and diagram
    P_inverter = 254.97  # Active power at the inverter in MW
    voltage_drop_percent = 2.5  # Voltage drop at Bus 6 in percent
    h5_distortion_percent = 10  # 5th harmonic distortion in percent
    h3_distortion_percent = 5   # 3rd harmonic distortion in percent

    # --- Step 1: Calculate reactive power compensation for the fault's impact ---
    Q_comp_fault = (voltage_drop_percent / 100) * P_inverter

    # --- Step 2: Calculate reactive power compensation for harmonic distortion ---
    total_harmonic_percent = h5_distortion_percent + h3_distortion_percent
    Q_comp_harmonics = (total_harmonic_percent / 100) * P_inverter

    # --- Step 3: Calculate the total reactive power compensation ---
    Q_total_compensation = Q_comp_fault + Q_comp_harmonics

    # --- Step 4: Print the detailed calculation ---
    print("Calculation of Total Reactive Power Compensation")
    print("="*50)

    print("Part 1: Compensation for Fault Impact (Voltage Drop)")
    print(f"The compensation is calculated as {voltage_drop_percent}% of the inverter active power.")
    print(f"Q_fault = ({voltage_drop_percent} / 100) * {P_inverter}")
    print(f"Q_fault = {Q_comp_fault:.2f} MVAr\n")

    print("Part 2: Compensation for Harmonic Distortion")
    print(f"The total harmonic distortion percentage is {h5_distortion_percent}% (5th) + {h3_distortion_percent}% (3rd) = {total_harmonic_percent}%.")
    print(f"The compensation is calculated as {total_harmonic_percent}% of the inverter active power.")
    print(f"Q_harmonics = ({total_harmonic_percent} / 100) * {P_inverter}")
    print(f"Q_harmonics = {Q_comp_harmonics:.2f} MVAr\n")

    print("Part 3: Total Reactive Power Compensation")
    print("The total compensation is the sum of the two parts.")
    print(f"Q_total = Q_fault + Q_harmonics")
    print(f"Q_total = {Q_comp_fault:.2f} MVAr + {Q_comp_harmonics:.2f} MVAr")
    print(f"Q_total = {Q_total_compensation:.2f} MVAr")
    print("="*50)
    
    # Print the final answer in the specified format
    print(f"\n<<<The final answer is>>>\n<<<{Q_total_compensation:.1f}>>>")

# Execute the calculation
calculate_reactive_power_compensation()