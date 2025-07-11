import math

def calculate_reactive_power_compensation():
    """
    Calculates the total reactive power compensation required for an HVDC system.

    The calculation considers two main factors:
    1. Compensation for a voltage drop caused by a fault.
    2. Compensation for harmonic distortion from the converter.
    """

    # --- Given Parameters ---
    # Assuming AC Line-to-Line voltage is the same as the HVDC voltage
    V_ac_ll = 500e3  # V
    # System inductance used to find short-circuit level
    L = 0.1         # H
    # Assumed standard system frequency
    f = 50          # Hz
    # Voltage drop to be compensated
    voltage_drop_ratio = 2.5 / 100
    # Power measurements at the inverter
    P_inverter = 254.97e6 # W (This is the fundamental active power, P1)
    Q_net_inverter = 14.58e6 # Var (This is the net fundamental reactive power, Q1)
    # Harmonic current distortion ratios relative to the fundamental
    h5_distortion_ratio = 10 / 100
    h3_distortion_ratio = 5 / 100

    # --- Step 1: Calculate Reactive Power Compensation for Voltage Drop ---
    print("--- Step 1: Calculate Reactive Power Compensation for Voltage Drop ---")

    # Calculate the AC system reactance (X)
    omega = 2 * math.pi * f
    X_reactance = omega * L
    print(f"System reactance X = 2 * pi * f * L = 2 * 3.14159 * {f} Hz * {L} H = {X_reactance:.2f} Ohms")

    # Calculate the short-circuit apparent power (S_sc)
    S_sc = V_ac_ll**2 / X_reactance
    S_sc_MVA = S_sc / 1e6
    print(f"Short-circuit power S_sc = V^2 / X = ({V_ac_ll/1e3:.0f} kV)^2 / {X_reactance:.2f} Ohms = {S_sc_MVA:.2f} MVA")

    # Calculate the reactive power compensation needed to restore the voltage
    Q_comp_fault = S_sc * voltage_drop_ratio
    Q_comp_fault_MVar = Q_comp_fault / 1e6
    print(f"Reactive power for fault compensation Q_fault = S_sc * (Delta_V/V) = {S_sc_MVA:.2f} MVA * {voltage_drop_ratio*100:.1f}% = {Q_comp_fault_MVar:.2f} MVar")
    print("\n" + "="*60 + "\n")

    # --- Step 2: Calculate Reactive Power Compensation for Harmonics ---
    print("--- Step 2: Calculate Reactive Power Compensation for Harmonics ---")
    print("(Calculated as the Distortion Reactive Power, D)")

    # The given power values are for the fundamental frequency
    P1 = P_inverter
    Q1 = Q_net_inverter

    # Calculate fundamental apparent power (S1)
    S1 = math.sqrt(P1**2 + Q1**2)
    S1_MVA = S1 / 1e6
    print(f"Fundamental apparent power S1 = sqrt(P1^2 + Q1^2) = sqrt(({P1/1e6:.2f} MW)^2 + ({Q1/1e6:.2f} MVar)^2) = {S1_MVA:.2f} MVA")

    # Calculate total harmonic distortion of current (THD_i)
    THD_i = math.sqrt(h3_distortion_ratio**2 + h5_distortion_ratio**2)
    print(f"Current THD_i = sqrt((I3/I1)^2 + (I5/I1)^2) = sqrt({h3_distortion_ratio:.2f}^2 + {h5_distortion_ratio:.2f}^2) = {THD_i:.4f}")

    # Calculate distortion reactive power (D)
    Q_comp_harmonics = THD_i * S1
    Q_comp_harmonics_MVar = Q_comp_harmonics / 1e6
    print(f"Reactive power for harmonic compensation Q_harmonics = S1 * THD_i = {S1_MVA:.2f} MVA * {THD_i:.4f} = {Q_comp_harmonics_MVar:.2f} MVar")
    print("\n" + "="*60 + "\n")

    # --- Step 3: Calculate Total Reactive Power Compensation ---
    print("--- Step 3: Calculate Total Reactive Power Compensation ---")

    # Sum the two components
    Q_total_comp_MVar = Q_comp_fault_MVar + Q_comp_harmonics_MVar
    print("Total reactive power compensation required is the sum of the compensation for the fault and for the harmonics.")
    print(f"Q_total = Q_fault + Q_harmonics")
    print(f"Q_total = {Q_comp_fault_MVar:.2f} MVar + {Q_comp_harmonics_MVar:.2f} MVar = {Q_total_comp_MVar:.2f} MVar")
    return Q_total_comp_MVar

if __name__ == '__main__':
    total_compensation = calculate_reactive_power_compensation()
    # The final numerical answer is requested in a specific format
    # print(f"\n<<<Total Compensation: {total_compensation:.2f} MVar>>>")
    print(f'<<<{total_compensation:.2f}>>>')