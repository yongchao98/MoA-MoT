import numpy as np

def calculate_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """
    # Step 1: Parse Parameters
    V_rf_fundamental = 1.0  # V (peak) from the table
    f0 = 915e6  # Hz (fundamental frequency) from the table (ω = 2πf) and parasitic model
    R0 = 50.0  # Ohms, base parasitic resistance
    C_parasitic = 2e-15  # Farads (2 fF)
    
    # Step 3: Assume the rectifier's input impedance is purely resistive for simplification.
    R_load = 50.0  # Ohms

    # Step 2: Model Harmonic Voltages
    harmonics = [1, 3, 5, 7]
    voltages = [
        V_rf_fundamental,
        V_rf_fundamental * 0.9,          # 10% drop from the 1st harmonic
        V_rf_fundamental * 0.9**2,       # 10% drop from the 3rd harmonic
        V_rf_fundamental * 0.9**3        # 10% drop from the 5th harmonic
    ]

    total_power_in = 0
    total_power_lost = 0
    total_power_delivered = 0

    power_delivered_list = []
    power_in_list = []

    print("Calculating system efficiency based on the provided model...")
    print("-" * 50)
    print("Assumptions:")
    print("1. Harmonic voltages: V_n = 0.9 * V_{n-2}, starting from V_1 = 1.0V.")
    print(f"2. Rectifier input impedance is modeled as a simple R_load = {R_load} Ohm resistor.")
    print("-" * 50)

    # Step 4: Iterate and Calculate Power per Harmonic
    for i, n in enumerate(harmonics):
        v_peak = voltages[i]
        f = n * f0
        omega = 2 * np.pi * f

        # Calculate frequency-dependent parasitic resistance
        r_parasitic = R0 * (n)**2

        # Calculate parasitic capacitive reactance
        x_c_parasitic = -1 / (omega * C_parasitic)

        # Total series impedance: Z_total = R_parasitic + R_load + j*X_c_parasitic
        z_real = r_parasitic + R_load
        z_imag = x_c_parasitic
        z_mag = np.sqrt(z_real**2 + z_imag**2)

        # Peak current for the harmonic
        i_peak = v_peak / z_mag

        # Average power calculation (P = 0.5 * I_peak^2 * R)
        power_lost = 0.5 * i_peak**2 * r_parasitic
        power_delivered = 0.5 * i_peak**2 * R_load
        power_in = power_lost + power_delivered

        # Store for final equation
        power_delivered_list.append(power_delivered)
        power_in_list.append(power_in)

        # Accumulate totals
        total_power_in += power_in
        total_power_lost += power_lost
        total_power_delivered += power_delivered
        
        print(f"For Harmonic n={n}:")
        print(f"  - Peak Voltage (V_peak): {v_peak:.3f} V")
        print(f"  - Frequency: {f/1e6:.1f} MHz")
        print(f"  - Parasitic Resistance (R_p): {r_parasitic:.1f} Ohm")
        print(f"  - Parasitic Reactance (X_c): {x_c_parasitic:.1f} Ohm")
        print(f"  - Total Impedance (|Z|): {z_mag:.1f} Ohm")
        print(f"  - Peak Current (I_peak): {i_peak*1e6:.3f} uA")
        print(f"  - Power Lost in R_p: {power_lost*1e6:.3f} uW")
        print(f"  - Power Delivered to Load: {power_delivered*1e6:.3f} uW")
        print("-" * 20)

    # Step 5: Calculate Overall Efficiency
    efficiency = total_power_delivered / total_power_in

    print("Summary of Power:")
    print(f"Total Power Input from Source = {total_power_in*1e6:.3f} uW")
    print(f"Total Power Lost in Parasitics = {total_power_lost*1e6:.3f} uW")
    print(f"Total Power Delivered to Load = {total_power_delivered*1e6:.3f} uW")
    print("\n" + "="*50)
    print("Final Efficiency Calculation")
    print("="*50)
    print("Efficiency = Total Power Delivered / Total Power Input")
    
    p_del_str = " + ".join([f"{p:.3e}" for p in power_delivered_list])
    p_in_str = " + ".join([f"{p:.3e}" for p in power_in_list])
    
    print(f"Efficiency = ({p_del_str}) W / ({p_in_str}) W")
    
    p_del_sum = sum(power_delivered_list)
    p_in_sum = sum(power_in_list)
    
    print(f"Efficiency = {p_del_sum:.3e} W / {p_in_sum:.3e} W")
    print(f"Overall System Efficiency = {efficiency*100:.2f}%")


calculate_efficiency()
<<<2.84>>>