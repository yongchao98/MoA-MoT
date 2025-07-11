import math

def calculate_reactive_power_compensation():
    """
    Calculates the total reactive power compensation required for an HVDC system.
    The calculation includes compensation for harmonic distortion and voltage stability after a fault.
    """
    # Step 1: Define the given parameters from the problem description.
    # Active power at the inverter bus (Bus 6), from the diagram.
    P_inverter = 254.97  # in MW
    # Voltage drop at Bus 6 after the fault.
    voltage_drop_percent = 2.5
    voltage_drop_ratio = voltage_drop_percent / 100
    # Harmonic distortion percentages.
    harmonic_dist_5th_percent = 10
    harmonic_dist_3rd_percent = 5
    harmonic_dist_5th_ratio = harmonic_dist_5th_percent / 100
    harmonic_dist_3rd_ratio = harmonic_dist_3rd_percent / 100
    # Constant for voltage stability reactive power estimation
    K_volt_stability = 2

    # Step 2: Calculate the reactive power required for harmonic filtering (Q_harm).
    # This is modeled as a percentage of the active power.
    total_harmonic_ratio = harmonic_dist_5th_ratio + harmonic_dist_3rd_ratio
    Q_harm = total_harmonic_ratio * P_inverter
    
    # Step 3: Calculate the reactive power required to restore voltage stability (Q_volt).
    # This is estimated based on the voltage drop and active power.
    Q_volt = K_volt_stability * voltage_drop_ratio * P_inverter

    # Step 4: Calculate the total reactive power compensation.
    Q_total = Q_harm + Q_volt
    
    # Print the explanation and the final calculation step-by-step.
    print("This script calculates the total reactive power compensation required for the HVDC system.")
    print("-" * 50)
    
    print("1. Calculate reactive power for harmonic filtering (Q_harm):")
    print(f"Q_harm = (Distortion_5th + Distortion_3rd) * P_inverter")
    print(f"Q_harm = ({harmonic_dist_5th_ratio} + {harmonic_dist_3rd_ratio}) * {P_inverter} MW")
    print(f"Q_harm = {total_harmonic_ratio:.2f} * {P_inverter} MW = {Q_harm:.2f} MVar")
    print("")

    print("2. Calculate reactive power for voltage stability (Q_volt):")
    print(f"Q_volt = k * (Voltage_Drop_Ratio) * P_inverter")
    print(f"Q_volt = {K_volt_stability} * {voltage_drop_ratio} * {P_inverter} MW")
    print(f"Q_volt = {K_volt_stability * voltage_drop_ratio:.2f} * {P_inverter} MW = {Q_volt:.2f} MVar")
    print("")

    print("3. Calculate total reactive power compensation (Q_total):")
    print(f"Q_total = Q_harm + Q_volt")
    print(f"Q_total = {Q_harm:.2f} MVar + {Q_volt:.2f} MVar")
    print(f"The total reactive power compensation required is {Q_total:.2f} MVar.")
    print("-" * 50)

    # Returning the final numerical answer for the platform
    return Q_total

if __name__ == "__main__":
    result = calculate_reactive_power_compensation()
    # The final answer in the required format
    # print(f"<<<{result:.2f}>>>")
    
# Call the function to execute the calculation and print the results.
calculate_reactive_power_compensation()