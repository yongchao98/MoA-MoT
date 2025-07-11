import math

def calculate_reactive_power_compensation():
    """
    Calculates the total reactive power compensation required for the HVDC system.
    """

    # --- Step 0: Define Given Parameters from the problem ---
    # Active power transmitted by the inverter at bus 6
    P_inv = 254.97  # MW
    # Voltage drop at bus 6 after the fault
    voltage_drop_pct = 2.5  # %
    # Harmonic distortion percentages
    h3_distortion_pct = 5.0 # % due to third harmonic
    h5_distortion_pct = 10.0 # % due to fifth harmonic

    print("Step-by-Step Calculation for Total Reactive Power Compensation:\n")

    # --- Step 1: Calculate the fundamental reactive power demand (Q_fund) ---
    # A common rule of thumb is that the reactive power consumed by an LCC converter
    # is approximately 50% of its active power transmission.
    q_fund_ratio = 0.5
    Q_fund = P_inv * q_fund_ratio
    print(f"1. Fundamental Reactive Power Demand (Q_fund):")
    print(f"   Q_fund = P_inverter * 50%")
    print(f"   Q_fund = {P_inv:.2f} MW * {q_fund_ratio:.2f} = {Q_fund:.2f} MVar\n")

    # --- Step 2: Calculate the harmonic reactive power demand (Q_harm) ---
    # This is calculated based on the harmonic currents.
    h3_rel = h3_distortion_pct / 100
    h5_rel = h5_distortion_pct / 100
    
    # Calculate the Total Harmonic Distortion of the current (THD_I)
    thd_i = math.sqrt(h3_rel**2 + h5_rel**2)
    
    # Calculate the fundamental apparent power (S_fund)
    S_fund = math.sqrt(P_inv**2 + Q_fund**2)
    
    # The harmonic reactive power (or distortion power) is approximated by S_fund * THD_I
    Q_harm = S_fund * thd_i
    print(f"2. Harmonic Reactive Power Demand (Q_harm):")
    print(f"   First, calculate current THD_I = sqrt({h3_rel:.2f}² + {h5_rel:.2f}²) = {thd_i:.4f}")
    print(f"   Then, find fundamental apparent power S_fund = sqrt({P_inv:.2f}² + {Q_fund:.2f}²) = {S_fund:.2f} MVA")
    print(f"   Q_harm = S_fund * THD_I = {S_fund:.2f} MVA * {thd_i:.4f} = {Q_harm:.2f} MVar\n")

    # --- Step 3: Calculate the reactive power deficit due to the voltage sag (Q_sag) ---
    # The voltage drop reduces the output of existing capacitor banks (Q is proportional to V^2).
    # This lost reactive power must be compensated for.
    voltage_factor = 1 - (voltage_drop_pct / 100)
    
    # The deficit is the difference in the capacitor's output before and after the sag.
    # Q_sag = Q_pre_fault * (1 - V_new^2/V_old^2)
    # Assume the pre-fault compensation matched the fundamental demand, so Q_pre_fault = Q_fund.
    Q_sag = Q_fund * (1 - voltage_factor**2)
    print(f"3. Reactive Power Deficit from Voltage Sag (Q_sag):")
    print(f"   Due to the {voltage_drop_pct}% voltage drop, existing compensation is less effective.")
    print(f"   Q_sag = Q_fund * (1 - (1 - {voltage_drop_pct/100:.3f})²)")
    print(f"   Q_sag = {Q_fund:.2f} MVar * (1 - {voltage_factor:.3f}²) = {Q_sag:.2f} MVar\n")

    # --- Step 4: Calculate the total required reactive power compensation (Q_total) ---
    Q_total = Q_fund + Q_harm + Q_sag
    print("4. Total Reactive Power Compensation Required (Q_total):")
    print("   The total compensation is the sum of the fundamental, harmonic, and sag-related components.")
    print(f"   Q_total = Q_fund + Q_harm + Q_sag")
    print(f"   Q_total = {Q_fund:.2f} MVar + {Q_harm:.2f} MVar + {Q_sag:.2f} MVar = {Q_total:.2f} MVar")


if __name__ == '__main__':
    calculate_reactive_power_compensation()