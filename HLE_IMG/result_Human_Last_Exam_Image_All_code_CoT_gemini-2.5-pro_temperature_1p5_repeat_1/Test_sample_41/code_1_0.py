import math

def calculate_reactive_power_compensation():
    """
    Calculates the total reactive power compensation required for an HVDC system.

    The calculation considers two factors:
    1. Reactive power needed to compensate for harmonic distortion.
    2. Reactive power needed to restore voltage stability after a fault.
    """

    # --- Given Parameters ---
    # Real power at the inverter (from the diagram)
    P_inverter = 254.97  # in MW

    # Harmonic distortion percentages
    hd5_percent = 10.0
    hd3_percent = 5.0

    # Voltage drop percentage after fault
    voltage_drop_percent = 2.5

    # --- Calculations ---

    # 1. Calculate compensation for harmonic distortion
    # Total percentage for harmonic compensation
    total_hd_percent = hd5_percent + hd3_percent
    # Reactive power for harmonics (Q_harmonics)
    Q_harmonics = (total_hd_percent / 100.0) * P_inverter

    # 2. Calculate compensation for voltage drop
    # Reactive power for voltage support (Q_voltage)
    Q_voltage_support = (voltage_drop_percent / 100.0) * P_inverter

    # 3. Calculate total reactive power compensation
    Q_total = Q_harmonics + Q_voltage_support

    # --- Output the results ---
    print("This script calculates the total reactive power compensation required.")
    print("-" * 60)
    print(f"Base Real Power (P_inverter): {P_inverter} MW")
    print(f"Harmonic Distortion: {hd5_percent}% (5th) + {hd3_percent}% (3rd) = {total_hd_percent}%")
    print(f"Post-Fault Voltage Drop: {voltage_drop_percent}%")
    print("-" * 60)
    
    print("Final Equation:")
    # The final equation requires printing each number
    print(f"Total Compensation = (({hd5_percent}/100 * {P_inverter}) + ({voltage_drop_percent}/100 * {P_inverter})) MVar")
    print(f"Total Compensation = {Q_harmonics:.4f} MVar + {Q_voltage_support:.4f} MVar")
    print(f"Total Compensation = {Q_total:.4f} MVar")

# Execute the function
calculate_reactive_power_compensation()

# For direct answer extraction
Q_harmonics = ((10.0 + 5.0) / 100.0) * 254.97
Q_voltage_support = (2.5 / 100.0) * 254.97
Q_total = Q_harmonics + Q_voltage_support
print(f"\n<<<{Q_total:.2f}>>>")