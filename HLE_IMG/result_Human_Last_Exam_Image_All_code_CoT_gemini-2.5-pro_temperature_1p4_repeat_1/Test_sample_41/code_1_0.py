import math

def calculate_reactive_power_compensation():
    """
    Calculates the total reactive power compensation required for an HVDC system.

    The calculation is based on two components:
    1. Reactive power to restore voltage after a fault.
    2. Reactive power to compensate for harmonic distortion from the converter.
    """

    # --- Given Parameters ---
    P_inv_MW = 254.97  # Active power at the inverter in MW
    V_hvdc_kV = 500.0   # HVDC line voltage in kV
    L_H = 0.1         # System inductance in Henrys
    voltage_drop_pct = 2.5 # Voltage drop percentage at Bus 6
    hd5_pct = 10.0      # 5th harmonic distortion percentage
    hd3_pct = 5.0       # 3rd harmonic distortion percentage
    f_Hz = 50.0         # Assumed system frequency in Hz

    # --- Step 1: Calculate Reactive Power for Voltage Support (Q_V) ---

    # Assume nominal AC voltage is equal to the HVDC voltage
    V_nom_V = V_hvdc_kV * 1000
    
    # Calculate the final voltage after the drop
    V_final_V = V_nom_V * (1 - voltage_drop_pct / 100.0)

    # Calculate the equivalent system reactance
    X_th_Ohm = 2 * math.pi * f_Hz * L_H

    # Calculate the reactive power required to restore voltage
    # Q_V = (V_nom^2 - V_final^2) / X_th
    Q_V_Var = (V_nom_V**2 - V_final_V**2) / X_th_Ohm
    Q_V_MVar = Q_V_Var / 1e6

    print("--- Part 1: Compensation for Voltage Drop ---")
    print(f"Nominal Voltage (V_nom): {V_nom_V/1000:.2f} kV")
    print(f"Voltage after drop (V_final): {V_final_V/1000:.2f} kV")
    print(f"System Reactance (X_th): {X_th_Ohm:.2f} Ohms")
    print(f"Reactive Power for Voltage Support (Q_V): ({V_nom_V**2:.2e} - {V_final_V**2:.2e}) / {X_th_Ohm:.2f} = {Q_V_MVar:.2f} MVar\n")


    # --- Step 2: Calculate Reactive Power for Harmonic Compensation (Q_H) ---

    # Calculate reactive power compensation based on active power and harmonic distortion percentages
    Q_H5_MVar = P_inv_MW * (hd5_pct / 100.0)
    Q_H3_MVar = P_inv_MW * (hd3_pct / 100.0)
    Q_H_MVar = Q_H5_MVar + Q_H3_MVar

    print("--- Part 2: Compensation for Harmonic Distortion ---")
    print(f"Active Power at Inverter (P_inv): {P_inv_MW:.2f} MW")
    print(f"Compensation for 5th Harmonic ({hd5_pct}% of P): {Q_H5_MVar:.2f} MVar")
    print(f"Compensation for 3rd Harmonic ({hd3_pct}% of P): {Q_H3_MVar:.2f} MVar")
    print(f"Total Reactive Power for Harmonics (Q_H): {Q_H5_MVar:.2f} + {Q_H3_MVar:.2f} = {Q_H_MVar:.2f} MVar\n")

    # --- Step 3: Calculate Total Reactive Power Compensation ---
    Q_total_MVar = Q_V_MVar + Q_H_MVar
    
    print("--- Part 3: Total Compensation ---")
    print("Total Reactive Power Compensation (Q_total) = Q_V + Q_H")
    print(f"Q_total = {Q_V_MVar:.2f} MVar + {Q_H_MVar:.2f} MVar")
    print(f"Q_total = {Q_total_MVar:.2f} MVar")
    
    return Q_total_MVar

if __name__ == '__main__':
    total_compensation = calculate_reactive_power_compensation()
    # The final numerical answer is wrapped in <<<>>> as requested.
    print(f"\nFinal calculated total reactive power compensation is {total_compensation:.2f} MVar.")
    # The final answer format requirement.
    final_answer = round(total_compensation, 2)
    # The final answer must be a single value, as per the example. I will choose the most final calculation result.
    print(f'<<<EOD>>>') # Not outputting the value directly but printing the final equation is the spirit.
    # The prompt actually said "directly return the answer with the format <<<answer content>>> at the end of your response". This implies one value. Let's do it.
    print(f'<<<{final_answer}>>>')
