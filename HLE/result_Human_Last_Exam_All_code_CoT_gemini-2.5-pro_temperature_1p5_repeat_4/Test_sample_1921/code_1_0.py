import math

def main():
    """
    Calculates the Signal-to-Noise Ratio (SNR) for a microwave link.
    """
    # --- Input Parameters ---
    # Transmitter (Tx)
    tx_power_dbm = 30.0  # dBm
    tx_antenna_gain_db = 20.0  # dB
    tx_antenna_loss_db = 1.0  # dB
    tx_filter_loss_db = 1.0  # dB
    tx_cable_loss_db = 1.0  # dB

    # Path
    frequency_ghz = 24.0  # GHz
    distance_km = 10.0  # km

    # Receiver (Rx)
    rx_antenna_gain_db = 1.0  # dB
    rx_antenna_loss_db = 0.5  # dB
    rx_input_filter_loss_db = 1.0 # dB
    lna_gain_db = 36.0 # dB
    lna_nf_db = 2.0 # dB
    mixer_loss_db = 9.0 # dB (conversion loss)
    mixer_nf_db = 9.0 # dB (NF of passive mixer equals its loss)
    if_filter_loss_db = 1.0 # dB

    # System
    bandwidth_hz = 100e3  # 100 kHz
    temperature_k = 300.0  # Kelvin
    k_boltzmann = 1.38e-23  # J/K

    # --- 1. Calculate EIRP (Effective Isotropic Radiated Power) ---
    total_tx_loss_db = tx_antenna_loss_db + tx_filter_loss_db + tx_cable_loss_db
    eirp_dbm = tx_power_dbm + tx_antenna_gain_db - total_tx_loss_db
    print(f"Step 1: Calculate Transmitter EIRP")
    print(f"EIRP (dBm) = Tx Power (dBm) + Tx Ant Gain (dB) - Tx Losses (dB)")
    print(f"EIRP (dBm) = {tx_power_dbm} + {tx_antenna_gain_db} - ({tx_antenna_loss_db} + {tx_filter_loss_db} + {tx_cable_loss_db}) = {eirp_dbm:.2f} dBm\n")

    # --- 2. Calculate Free Space Path Loss (FSPL) ---
    # Using formula: FSPL(dB) = 32.45 + 20*log10(f_MHz) + 20*log10(d_km)
    frequency_mhz = frequency_ghz * 1000
    fspl_db = 32.45 + 20 * math.log10(frequency_mhz) + 20 * math.log10(distance_km)
    print(f"Step 2: Calculate Free Space Path Loss (FSPL)")
    print(f"FSPL (dB) = 32.45 + 20*log10({frequency_mhz:.0f} MHz) + 20*log10({distance_km:.0f} km) = {fspl_db:.2f} dB\n")

    # --- 3. Calculate Received Signal Power (S) ---
    # This is the power at the input of the receiver chain (after Rx antenna gain).
    signal_power_dbm = eirp_dbm - fspl_db + rx_antenna_gain_db
    print(f"Step 3: Calculate Received Signal Power (S)")
    print(f"S (dBm) = EIRP (dBm) - FSPL (dB) + Rx Ant Gain (dB)")
    print(f"S (dBm) = {eirp_dbm:.2f} - {fspl_db:.2f} + {rx_antenna_gain_db} = {signal_power_dbm:.2f} dBm\n")

    # --- 4. Calculate Total System Noise Power (N) ---
    print(f"Step 4: Calculate Total System Noise Power (N)")
    # Step 4a: Calculate Cascaded System Noise Figure (NF) using Friis formula
    # F_sys = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
    # We need to convert gains and NFs to linear scale.
    # Note: For passive components at ambient temp, Loss (dB) = NF (dB)
    g1_lin = 10**(-rx_antenna_loss_db / 10)
    f1_lin = 10**(rx_antenna_loss_db / 10)
    
    g2_lin = 10**(-rx_input_filter_loss_db / 10)
    f2_lin = 10**(rx_input_filter_loss_db / 10)

    g3_lin = 10**(lna_gain_db / 10)
    f3_lin = 10**(lna_nf_db / 10)

    g4_lin = 10**(-mixer_loss_db / 10)
    f4_lin = 10**(mixer_nf_db / 10)

    # Subsequent stages will have negligible contribution due to high LNA gain
    f_sys_lin = f1_lin + (f2_lin - 1) / g1_lin + (f3_lin - 1) / (g1_lin * g2_lin) + \
                (f4_lin - 1) / (g1_lin * g2_lin * g3_lin)
    nf_sys_db = 10 * math.log10(f_sys_lin)
    print(f" - System Noise Figure (NF) calculated using Friis formula = {nf_sys_db:.2f} dB")
    
    # Step 4b: Calculate Thermal Noise Power
    thermal_noise_power_watts = k_boltzmann * temperature_k * bandwidth_hz
    thermal_noise_power_dbm = 10 * math.log10(thermal_noise_power_watts / 0.001)
    print(f" - Thermal Noise Power in {bandwidth_hz/1000:.0f} kHz BW = {thermal_noise_power_dbm:.2f} dBm")
    
    # Step 4c: Calculate Total Noise Power
    total_noise_power_dbm = thermal_noise_power_dbm + nf_sys_db
    print(f" - Total Noise Power (N) = Thermal Noise (dBm) + System NF (dB)")
    print(f"   N (dBm) = {thermal_noise_power_dbm:.2f} + {nf_sys_db:.2f} = {total_noise_power_dbm:.2f} dBm\n")

    # --- 5. Calculate Final SNR ---
    snr_db = signal_power_dbm - total_noise_power_dbm
    print(f"--- FINAL CALCULATION ---")
    print(f"SNR (dB) = Signal Power (dBm) - Total Noise Power (dBm)")
    print(f"SNR (dB) = {signal_power_dbm:.2f} - ({total_noise_power_dbm:.2f})")
    print(f"Final SNR = {snr_db:.2f} dB")

if __name__ == '__main__':
    main()