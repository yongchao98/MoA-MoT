import math

def calculate_reactive_power_compensation():
    """
    Calculates the total reactive power compensation required for an HVDC system.

    The calculation is based on:
    1. Compensation for a voltage drop after a fault.
    2. Compensation for harmonic distortions.
    """

    # --- Given values from the problem description and diagram ---
    # Power at inverter (Bus 6)
    P_inv = 254.97  # MW
    Q_inv_net = 14.58  # MVar (net reactive power consumed at bus 6)

    # System conditions
    voltage_drop_percentage = 2.5  # %
    harmonic_dist_3rd = 5.0      # %
    harmonic_dist_5th = 10.0     # %

    # --- Calculations ---

    print("Step 1: Calculate Reactive Power Compensation for Voltage Drop")
    # Convert voltage drop percentage to a decimal
    V_drop = voltage_drop_percentage / 100.0

    # Calculate the apparent power (S) at the inverter bus
    S_inv = math.sqrt(P_inv**2 + Q_inv_net**2)
    print(f"The apparent power at the inverter bus (S_inv) is calculated as: sqrt({P_inv}^2 + {Q_inv_net}^2) = {S_inv:.2f} MVA")

    # Calculate the compensation required to restore the voltage.
    # The formula used is Q_voltage = S * (1 - V_final^2), where V_final is the post-fault voltage in per unit.
    Q_voltage = S_inv * (1 - (1 - V_drop)**2)
    print(f"The reactive power needed for voltage support (Q_voltage) is: {S_inv:.2f} MVA * (1 - (1 - {V_drop})^2) = {Q_voltage:.2f} MVar\n")

    print("Step 2: Calculate Reactive Power Compensation for Harmonics")
    # Convert harmonic distortion percentages to decimals
    D3 = harmonic_dist_3rd / 100.0
    D5 = harmonic_dist_5th / 100.0

    # The required compensation is estimated as a percentage of the active power.
    Q_harmonic = (D3 + D5) * P_inv
    print(f"The reactive power for harmonic compensation (Q_harmonic) is estimated based on active power:")
    print(f"Q_harmonic = ({D3} + {D5}) * {P_inv:.2f} MW = {Q_harmonic:.2f} MVar\n")

    print("Step 3: Calculate Total Reactive Power Compensation")
    # The total compensation is the sum of the two parts.
    Q_total = Q_voltage + Q_harmonic

    # Print the final equation with all the calculated numbers
    print("The final equation for total reactive power compensation is:")
    print(f"Total Compensation = Reactive Power for Voltage Support + Reactive Power for Harmonics")
    print(f"Total Compensation = {Q_voltage:.2f} MVar + {Q_harmonic:.2f} MVar = {Q_total:.2f} MVar")

if __name__ == '__main__':
    calculate_reactive_power_compensation()