import math

def calculate_reactive_power_compensation():
    """
    Calculates the total reactive power compensation required for an HVDC system.

    This function considers reactive power demand from harmonic distortion and the
    additional reactive power needed to restore voltage after a fault.
    """

    # --- Given System Parameters ---
    P_inverter = 254.97  # Active power at the inverter (MW)
    HD_5 = 0.10          # 10% harmonic distortion due to 5th harmonics
    HD_3 = 0.05          # 5% harmonic distortion due to 3rd harmonics
    V_drop_percentage = 2.5 # Voltage drop at Bus 6 in percent

    # --- Assumptions ---
    # Assumption 1: Reactive power compensation for harmonics is a direct percentage
    # of the active power transferred.
    # Assumption 2: Capacitor banks for reactive power support are sized to provide
    # 50% of the active power rating of the converter.
    Q_cap_sizing_factor = 0.5

    print("--- Calculating Reactive Power Compensation ---")

    # Step 1: Calculate the reactive power compensation required for harmonics (Q_h)
    print("\nStep 1: Calculate Reactive Power for Harmonics (Q_h)")
    total_harmonic_contribution = HD_3 + HD_5
    Q_h = total_harmonic_contribution * P_inverter
    print(f"The reactive power needed to compensate for harmonics is calculated as a percentage of the active power.")
    print(f"Q_h = (Distortion_3rd + Distortion_5th) * P_inverter")
    print(f"Q_h = ({HD_3:.2f} + {HD_5:.2f}) * {P_inverter} = {Q_h:.4f} MVAr")

    # Step 2: Calculate reactive power compensation required to restore voltage (Q_v)
    print("\nStep 2: Calculate Reactive Power for Voltage Stability (Q_v)")
    # Estimate the initial size of the capacitor bank
    Q_cap_initial = Q_cap_sizing_factor * P_inverter
    print(f"First, we estimate the size of the existing capacitor bank. It is assumed to be 50% of the active power.")
    print(f"Initial Capacitor Power (Q_cap_initial) = {Q_cap_sizing_factor:.2f} * {P_inverter} MW = {Q_cap_initial:.4f} MVAr")

    # Calculate the compensation needed due to the voltage drop
    V_final_ratio = 1 - (V_drop_percentage / 100.0)
    Q_v = Q_cap_initial * (1 - V_final_ratio**2)
    print(f"A voltage drop of {V_drop_percentage}% reduces the capacitor's output, as Q is proportional to V^2.")
    print(f"Compensation for Voltage Drop (Q_v) = Q_cap_initial * (1 - (V_new/V_old)^2)")
    print(f"Q_v = {Q_cap_initial:.4f} MVAr * (1 - {V_final_ratio:.3f}^2) = {Q_v:.4f} MVAr")

    # Step 3: Calculate the total reactive power compensation required
    print("\nStep 3: Calculate Total Required Compensation (Q_total)")
    Q_total = Q_h + Q_v
    print(f"The total compensation is the sum of the harmonic and voltage stability components.")
    print(f"Q_total = Q_h + Q_v")
    print(f"Q_total = {Q_h:.4f} MVAr + {Q_v:.4f} MVAr = {Q_total:.4f} MVAr")

    print("\n--- Final Answer ---")
    print(f"The total reactive power compensation required from the system is {Q_total:.2f} MVAr.")

# Execute the calculation
if __name__ == '__main__':
    calculate_reactive_power_compensation()
    # The final numerical answer is derived from the calculation.
    # Q_h = (0.05 + 0.10) * 254.97 = 38.2455
    # Q_cap_initial = 0.5 * 254.97 = 127.485
    # Q_v = 127.485 * (1 - (1 - 0.025)^2) = 6.2936
    # Q_total = 38.2455 + 6.2936 = 44.5391
    # Rounded to two decimal places: 44.54
    # final_answer = 44.54
    # print(f"<<<{final_answer}>>>")