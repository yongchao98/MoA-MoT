import math

def solve_reactive_power_compensation():
    """
    Calculates the total reactive power compensation required for the HVDC system.
    """
    # --- Given Parameters ---
    V_dc = 500e3  # HVDC Voltage, V
    L = 0.1  # System Inductance, H
    P_inverter = 254.97  # Inverter real power, MW
    voltage_drop_percent = 2.5 / 100  # Voltage drop at Bus 6, dimensionless
    h5_distortion_percent = 10.0 / 100  # 5th harmonic distortion, dimensionless
    h3_distortion_percent = 5.0 / 100  # 3rd harmonic distortion, dimensionless

    # --- Assumptions ---
    # Assume the AC side voltage base is the same as the DC voltage, as no AC voltage is specified.
    V_base = V_dc  # V
    # Assume standard system frequency
    f = 50  # Hz

    print("Step 1: Calculate Reactive Power Compensation for Harmonics (Q_harmonics)")

    # Reactive power for 5th harmonic
    Q_h5 = h5_distortion_percent * P_inverter
    print(f"Reactive power for 5th harmonic (10% of P_inverter) = {h5_distortion_percent:.2f} * {P_inverter:.2f} MW = {Q_h5:.2f} MVar")

    # Reactive power for 3rd harmonic
    Q_h3 = h3_distortion_percent * P_inverter
    print(f"Reactive power for 3rd harmonic (5% of P_inverter) = {h3_distortion_percent:.2f} * {P_inverter:.2f} MW = {Q_h3:.2f} MVar")

    # Total harmonic reactive power
    Q_harmonics = Q_h5 + Q_h3
    print(f"Total harmonic reactive power Q_harmonics = {Q_h5:.2f} MVar + {Q_h3:.2f} MVar = {Q_harmonics:.2f} MVar")
    print("-" * 50)

    print("Step 2: Calculate Reactive Power Compensation for Voltage Drop (Q_fault)")

    # Calculate system Thevenin reactance from the given inductance
    omega = 2 * math.pi * f
    X_th = omega * L
    print(f"System reactance X_th = 2 * pi * f * L = 2 * pi * {f} Hz * {L} H = {X_th:.2f} Ohms")

    # Calculate the reactive power injection needed to raise the voltage by 2.5%
    # Formula: Q_comp ≈ (ΔV/V) * V_base^2 / X_th
    # where ΔV/V is the fractional voltage change required.
    # The result V_base^2 / X_th is in Var, so we divide by 1e6 to get MVar.
    Q_fault_MVar = (voltage_drop_percent * V_base**2) / X_th / 1e6
    print(f"Reactive power for fault Q_fault = (ΔV/V * V_base^2) / X_th")
    print(f"Q_fault = ({voltage_drop_percent:.3f} * ({V_base/1e3:.0f} kV)^2) / {X_th:.2f} Ohms = {Q_fault_MVar:.2f} MVar")
    print("-" * 50)


    print("Step 3: Calculate Total Reactive Power Compensation (Q_total)")
    # Total reactive power compensation
    Q_total = Q_fault_MVar + Q_harmonics
    print(f"Total compensation Q_total = Q_fault + Q_harmonics")
    print(f"Final Equation: Q_total = {Q_fault_MVar:.2f} MVar + {Q_harmonics:.2f} MVar")
    print(f"Total Reactive Power Required = {Q_total:.2f} MVar")

    return Q_total

if __name__ == '__main__':
    final_answer = solve_reactive_power_compensation()
    # The final answer is wrapped according to the instruction format.
    # print(f"\n<<<{final_answer:.2f}>>>")