import math

def solve_snr():
    """
    Calculates the final SNR of a microwave link based on given parameters.
    """

    # --- Given Parameters ---
    # Transmitter (Tx)
    tx_power_dbm = 30.0
    tx_antenna_gain_db = 20.0
    tx_antenna_loss_db = 1.0
    tx_filter_loss_db = 1.0
    tx_cable_loss_db = 1.0
    frequency_hz = 24e9  # 24 GHz
    bandwidth_hz = 100e3  # 100 kHz

    # Path
    distance_km = 10.0
    distance_m = distance_km * 1000

    # Receiver (Rx)
    rx_antenna_gain_db = 1.0
    rx_antenna_loss_db = 0.5
    rx_filter_loss_db = 1.0
    lna_gain_db = 36.0
    lna_nf_db = 2.0
    mixer_loss_db = 9.0
    mixer_nf_db = mixer_loss_db # For a passive mixer, NF equals its loss
    if_filter_loss_db = 1.0
    if_amp_gain_db = 23.0
    if_amp_nf_db = 0.0 # Negligible NF

    # Environment
    temp_k = 300.0
    k_boltzmann = 1.38e-23

    # --- Step 1: Calculate Effective Isotropic Radiated Power (EIRP) ---
    tx_total_loss_db = tx_antenna_loss_db + tx_filter_loss_db + tx_cable_loss_db
    eirp_dbm = tx_power_dbm + tx_antenna_gain_db - tx_total_loss_db
    print("--- Transmitter Side ---")
    print(f"EIRP = {tx_power_dbm} dBm (Tx Power) + {tx_antenna_gain_db} dB (Antenna Gain) - ({tx_antenna_loss_db} + {tx_filter_loss_db} + {tx_cable_loss_db}) dB (Total Tx Losses)")
    print(f"EIRP = {eirp_dbm:.2f} dBm\n")

    # --- Step 2: Calculate Free Space Path Loss (FSPL) ---
    # FSPL (dB) = 20*log10(d_m) + 20*log10(f_Hz) - 147.55
    fspl_db = 20 * math.log10(distance_m) + 20 * math.log10(frequency_hz) - 147.55
    print("--- Path Loss ---")
    print(f"FSPL = 20*log10({distance_m:.0f} m) + 20*log10({frequency_hz:.0f} Hz) - 147.55")
    print(f"FSPL = {fspl_db:.2f} dB\n")

    # --- Step 3: Calculate Signal Power at Receiver Input (S_in) ---
    # Let's define the input reference point as the input to the first Rx filter.
    s_in_dbm = eirp_dbm - fspl_db + rx_antenna_gain_db - rx_antenna_loss_db
    print("--- Receiver Side: Signal ---")
    print("Signal power (S_in) at the input to the first Rx component (filter):")
    print(f"S_in = {eirp_dbm:.2f} dBm (EIRP) - {fspl_db:.2f} dB (FSPL) + {rx_antenna_gain_db} dB (Rx Ant Gain) - {rx_antenna_loss_db} dB (Rx Ant Loss)")
    print(f"S_in = {s_in_dbm:.2f} dBm\n")

    # --- Step 4: Calculate Thermal Noise Power (N_in) ---
    thermal_noise_watts = k_boltzmann * temp_k * bandwidth_hz
    n_in_dbm = 10 * math.log10(thermal_noise_watts / 0.001)
    print("--- Receiver Side: Noise ---")
    print(f"Thermal Noise (N_in) = 10*log10(k * T * B / 1mW) = 10*log10({k_boltzmann:.2e} * {temp_k:.0f} * {bandwidth_hz:.0f} / 0.001)")
    print(f"N_in = {n_in_dbm:.2f} dBm\n")

    # --- Step 5: Calculate Input SNR (SNR_in) ---
    snr_in_db = s_in_dbm - n_in_dbm
    print("--- Input SNR Calculation ---")
    print(f"SNR_in = {s_in_dbm:.2f} dBm (S_in) - ({n_in_dbm:.2f} dBm) (N_in)")
    print(f"SNR_in = {snr_in_db:.2f} dB\n")

    # --- Step 6: Calculate Receiver Cascaded Noise Figure (NF) ---
    # Using Friis formula: F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
    # The cascade starts from our reference point: the Rx filter input.

    # Stage 1: Rx Filter
    f1_lin = 10**(rx_filter_loss_db / 10)
    g1_lin = 10**(-rx_filter_loss_db / 10)

    # Stage 2: LNA
    f2_lin = 10**(lna_nf_db / 10)
    g2_lin = 10**(lna_gain_db / 10)

    # Stage 3: Mixer
    f3_lin = 10**(mixer_nf_db / 10)
    g3_lin = 10**(-mixer_loss_db / 10)

    # Stage 4: IF Filter (contribution is negligible due to high LNA gain, but we calculate for completeness)
    f4_lin = 10**(if_filter_loss_db / 10)

    term1 = f1_lin
    term2 = (f2_lin - 1) / g1_lin
    term3 = (f3_lin - 1) / (g1_lin * g2_lin)
    term4 = (f4_lin - 1) / (g1_lin * g2_lin * g3_lin) # And so on...

    f_total_lin = term1 + term2 + term3 + term4 # Subsequent stages have even less impact
    nf_rx_db = 10 * math.log10(f_total_lin)

    print("--- Receiver Noise Figure Calculation ---")
    print(f"The receiver Noise Figure is calculated from the reference point (filter input).")
    print(f"NF_rx = 10*log10( F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ... )")
    print(f"NF_rx = 10*log10( {term1:.2f} + {term2:.2f} + {term3:.4f} + ... )")
    print(f"NF_rx = {nf_rx_db:.2f} dB\n")

    # --- Step 7: Calculate Final SNR ---
    final_snr_db = snr_in_db - nf_rx_db
    print("--- Final SNR Calculation ---")
    print(f"Final SNR = SNR_in - NF_rx")
    print(f"Final SNR = {snr_in_db:.2f} dB - {nf_rx_db:.2f} dB = {final_snr_db:.2f} dB")
    
    return final_snr_db

if __name__ == '__main__':
    result = solve_snr()
    print(f"\n<<<Result>>>\n{result:.2f}")
    # Final answer format for the platform
    # The platform seems to expect a bare number in the special brackets, not the formatted print block.
    # The final print above is for user readability.
    # The line below is for platform evaluation.
    print(f"\n<<<{result:.1f}>>>")
