import math

def calculate_reactive_power_compensation():
    """
    Calculates the total reactive power compensation required for an HVDC system.

    The calculation is based on two components:
    1. Compensation for harmonic distortion from the inverter.
    2. Compensation to restore voltage stability after a fault.
    """

    # --- Given Parameters and Assumptions ---
    P_inv = 254.97  # Active power at the inverter (MW)
    voltage_drop_percent = 2.5  # Voltage drop at Bus 6 (%)
    d3_percent = 5.0  # 3rd harmonic distortion (%)
    d5_percent = 10.0 # 5th harmonic distortion (%)
    L_sys = 0.1  # Assumed AC system Thevenin inductance (H)
    V_ac_LL = 500.0  # Assumed AC line-to-line voltage (kV)
    f = 50.0  # Assumed system frequency (Hz)
    cos_phi_1 = 0.9  # Assumed fundamental power factor of the inverter

    print("--- Problem Data and Assumptions ---")
    print(f"Inverter Active Power (P_inv): {P_inv} MW")
    print(f"Voltage Drop: {voltage_drop_percent}%")
    print(f"3rd Harmonic Distortion: {d3_percent}%")
    print(f"5th Harmonic Distortion: {d5_percent}%")
    print(f"Assumed System Inductance (L): {L_sys} H")
    print(f"Assumed AC Line Voltage (V_LL): {V_ac_LL} kV")
    print(f"Assumed System Frequency (f): {f} Hz")
    print(f"Assumed Inverter Power Factor (cos(phi)): {cos_phi_1}")
    print("-" * 40)

    # --- Step 1: Calculate Harmonic Reactive Power Compensation (Q_h) ---

    print("Step 1: Calculating Harmonic Compensation (Q_h)")
    # Convert distortion percentages to per-unit values
    d3 = d3_percent / 100.0
    d5 = d5_percent / 100.0

    # Calculate fundamental apparent power (S1)
    S1 = P_inv / cos_phi_1

    # Calculate the square of the Total Harmonic Current Distortion (THD_I)
    THD_I_sq = d3**2 + d5**2

    # Calculate Distortion Power (Q_h)
    Q_h = S1 * math.sqrt(THD_I_sq)
    print(f"  Fundamental Apparent Power (S1) = {P_inv:.2f} MW / {cos_phi_1} = {S1:.2f} MVA")
    print(f"  THD_I_sq = {d3:.2f}^2 + {d5:.2f}^2 = {THD_I_sq:.4f}")
    print(f"  Harmonic Compensation (Q_h) = {S1:.2f} MVA * sqrt({THD_I_sq:.4f}) = {Q_h:.2f} MVar")
    print("-" * 40)

    # --- Step 2: Calculate Fault Reactive Power Compensation (Q_fault) ---

    print("Step 2: Calculating Fault Compensation (Q_fault)")
    # Calculate the AC system reactance (X)
    X = 2 * math.pi * f * L_sys

    # Convert kV to V for calculation
    V_ac_LL_volts = V_ac_LL * 1000

    # Convert voltage drop percentage to per-unit
    delta_V_per_unit = voltage_drop_percent / 100.0

    # Calculate the required reactive power compensation
    Q_fault_var = delta_V_per_unit * (V_ac_LL_volts**2) / X
    Q_fault = Q_fault_var / 1e6  # Convert from Var to MVar

    print(f"  System Reactance (X) = 2 * pi * {f:.1f} Hz * {L_sys:.1f} H = {X:.2f} Ohms")
    print(f"  Fault Compensation (Q_fault) = {delta_V_per_unit:.3f} * ({V_ac_LL:.1f} kV)^2 / {X:.2f} Ohms = {Q_fault:.2f} MVar")
    print("-" * 40)

    # --- Step 3: Calculate Total Reactive Power Compensation ---

    print("Step 3: Calculating Total Compensation (Q_total)")
    Q_total = Q_h + Q_fault
    print("  The final equation is: Q_total = Q_h + Q_fault")
    print(f"  Substituting the calculated values:")
    print(f"  Q_total = {Q_h:.2f} MVar + {Q_fault:.2f} MVar")
    print(f"  Total Reactive Power Compensation = {Q_total:.2f} MVar")

    # Final answer in the requested format
    print(f"<<<{Q_total:.2f}>>>")

# Run the calculation
calculate_reactive_power_compensation()