import math

def solve_snr():
    """
    Calculates the Signal-to-Noise Ratio (SNR) for a microwave link based on given parameters.
    """
    # --- Given Parameters ---
    # Transmitter (Tx)
    tx_power_dbm = 30.0
    tx_ant_gain_db = 20.0
    tx_ant_loss_db = 1.0
    tx_filter_loss_db = 1.0
    tx_cable_loss_db = 1.0

    # Path
    frequency_hz = 24e9  # 24 GHz
    distance_m = 10e3   # 10 km

    # Receiver (Rx)
    rx_ant_gain_db = 1.0
    rx_ant_loss_db = 0.5
    rx_input_filter_loss_db = 1.0
    
    lna_gain_db = 36.0
    lna_nf_db = 2.0
    
    mixer_loss_db = 9.0 # Conversion Loss
    
    if_filter_loss_db = 1.0
    
    if_amp_gain_db = 23.0
    if_amp_nf_db = 0.0 # Negligible
    
    # System
    temperature_k = 300.0
    bandwidth_hz = 100e3 # 100 kHz

    # --- Constants ---
    BOLTZMANN_K = 1.38e-23  # J/K

    # --- Step 1: Calculate Effective Isotropic Radiated Power (EIRP) ---
    tx_total_loss_db = tx_ant_loss_db + tx_filter_loss_db + tx_cable_loss_db
    eirp_dbm = tx_power_dbm + tx_ant_gain_db - tx_total_loss_db
    print("--- Calculation Steps ---")
    print("1. Effective Isotropic Radiated Power (EIRP)")
    print(f"   EIRP = Tx Power + Tx Antenna Gain - Total Tx Losses")
    print(f"   EIRP = {tx_power_dbm} dBm + {tx_ant_gain_db} dB - ({tx_ant_loss_db} dB + {tx_filter_loss_db} dB + {tx_cable_loss_db} dB) = {eirp_dbm:.2f} dBm\n")

    # --- Step 2: Calculate Free Space Path Loss (FSPL) ---
    # Using the formula FSPL (dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
    # where d is in meters, f is in Hz, c is 3e8 m/s
    fspl_db = 20 * math.log10(distance_m) + 20 * math.log10(frequency_hz) - 147.55
    print("2. Free Space Path Loss (FSPL)")
    print(f"   FSPL = 20*log10(distance_m) + 20*log10(freq_hz) - 147.55")
    print(f"   FSPL = 20*log10({distance_m:.0f}) + 20*log10({frequency_hz:.0f}) - 147.55 = {fspl_db:.2f} dB\n")

    # --- Step 3: Calculate Signal Power at LNA Input ---
    rx_pre_lna_loss_db = rx_ant_loss_db + rx_input_filter_loss_db
    signal_at_lna_input_dbm = eirp_dbm - fspl_db + rx_ant_gain_db - rx_pre_lna_loss_db
    print("3. Signal Power at LNA Input (S)")
    print(f"   S = EIRP - FSPL + Rx Antenna Gain - Rx Pre-LNA Losses")
    print(f"   S = {eirp_dbm:.2f} dBm - {fspl_db:.2f} dB + {rx_ant_gain_db} dB - ({rx_ant_loss_db} dB + {rx_input_filter_loss_db} dB) = {signal_at_lna_input_dbm:.2f} dBm\n")

    # --- Step 4: Calculate Thermal Noise Power ---
    thermal_noise_watts = BOLTZMANN_K * temperature_k * bandwidth_hz
    thermal_noise_dbm = 10 * math.log10(thermal_noise_watts / 0.001)
    print("4. Thermal Noise Power (kTB)")
    print(f"   Noise_thermal = 10*log10(k * T * B * 1000)")
    print(f"   Noise_thermal = 10*log10({BOLTZMANN_K:.2e} * {temperature_k:.0f} * {bandwidth_hz:.0f} * 1000) = {thermal_noise_dbm:.2f} dBm\n")

    # --- Step 5: Calculate Receiver Noise Figure (at LNA input) using Friis Formula ---
    # Convert dB values to linear scale for the formula
    # F_sys = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
    g_lna = 10**(lna_gain_db / 10)
    f_lna = 10**(lna_nf_db / 10)
    
    # For a passive component like a mixer or filter, NF_dB = Loss_dB
    g_mixer = 10**(-mixer_loss_db / 10)
    f_mixer = 10**(mixer_loss_db / 10)

    g_if_filter = 10**(-if_filter_loss_db / 10)
    f_if_filter = 10**(if_filter_loss_db / 10)

    # Calculate cascaded noise figure from LNA onwards
    f_rx_linear = f_lna + (f_mixer - 1)/g_lna + (f_if_filter - 1)/(g_lna * g_mixer)
    # The contribution from stages after the IF filter is negligible due to the high preceding gain
    rx_nf_db = 10 * math.log10(f_rx_linear)
    
    print("5. Receiver Cascaded Noise Figure (NF_rx at LNA input)")
    print(f"   F_rx_linear = F_lna + (F_mixer - 1)/G_lna + ...")
    print(f"   F_lna = 10^({lna_nf_db}/10) = {f_lna:.2f}")
    print(f"   G_lna = 10^({lna_gain_db}/10) = {g_lna:.0f}")
    print(f"   F_mixer = 10^({mixer_loss_db}/10) = {f_mixer:.2f}")
    print(f"   NF_rx = 10*log10({f_lna:.2f} + ({f_mixer:.2f} - 1)/{g_lna:.0f} + ...) = {rx_nf_db:.2f} dB\n")

    # --- Step 6: Calculate Total Noise Power at LNA Input ---
    total_noise_dbm = thermal_noise_dbm + rx_nf_db
    print("6. Total Noise Power at LNA Input (N)")
    print(f"   N = Thermal Noise Power + Receiver Noise Figure")
    print(f"   N = {thermal_noise_dbm:.2f} dBm + {rx_nf_db:.2f} dB = {total_noise_dbm:.2f} dBm\n")

    # --- Step 7: Calculate Final SNR ---
    snr_db = signal_at_lna_input_dbm - total_noise_dbm
    print("7. Final Signal-to-Noise Ratio (SNR)")
    print(f"   SNR = Signal Power (S) - Total Noise Power (N)")
    print(f"   SNR = {signal_at_lna_input_dbm:.2f} dBm - ({total_noise_dbm:.2f} dBm) = {snr_db:.2f} dB")
    
    return snr_db

if __name__ == '__main__':
    final_snr = solve_snr()
    # Final answer format as requested by the user
    print(f"\n<<<SNR = {final_snr:.2f} dB>>>")
