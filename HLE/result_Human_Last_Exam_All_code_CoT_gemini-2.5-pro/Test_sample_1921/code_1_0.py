import math

def solve_snr():
    """
    Calculates the SNR of a microwave link based on given parameters.
    """

    # --- Given Parameters ---
    # Transmitter (Tx)
    tx_power_dbm = 30.0
    tx_ant_gain_db = 20.0
    tx_ant_loss_db = 1.0
    tx_filter_loss_db = 1.0
    tx_cable_loss_db = 1.0
    freq_ghz = 24.0
    bw_khz = 100.0

    # Path
    distance_km = 10.0

    # Receiver (Rx)
    rx_ant_gain_db = 1.0
    rx_ant_loss_db = 0.5
    rx_filter_loss_db = 1.0
    lna_gain_db = 36.0
    lna_nf_db = 2.0
    mixer_loss_db = 9.0
    if_filter_loss_db = 1.0
    if_amp_gain_db = 23.0
    if_amp_nf_db = 0.0   # Negligible NF

    # Environment
    temp_k = 300.0
    k_boltzmann = 1.38e-23  # J/K

    # --- Helper functions ---
    def db_to_linear(db_val):
        return 10**(db_val / 10.0)

    def linear_to_db(lin_val):
        return 10 * math.log10(lin_val)

    print("Calculating SNR Step-by-Step:\n")

    # --- Step 1: Calculate EIRP ---
    print("Step 1: Calculate Transmitter EIRP")
    total_tx_loss_db = tx_ant_loss_db + tx_filter_loss_db + tx_cable_loss_db
    eirp_dbm = tx_power_dbm + tx_ant_gain_db - total_tx_loss_db
    print(f"EIRP = Tx Power (dBm) + Tx Antenna Gain (dB) - Tx Losses (dB)")
    print(f"EIRP = {tx_power_dbm:.1f} + {tx_ant_gain_db:.1f} - ({tx_ant_loss_db:.1f} + {tx_filter_loss_db:.1f} + {tx_cable_loss_db:.1f}) = {eirp_dbm:.2f} dBm\n")

    # --- Step 2: Calculate FSPL ---
    print("Step 2: Calculate Free Space Path Loss (FSPL)")
    fspl_db = 20 * math.log10(distance_km) + 20 * math.log10(freq_ghz) + 92.45
    print(f"FSPL = 20*log10(d_km) + 20*log10(f_GHz) + 92.45")
    print(f"FSPL = 20*log10({distance_km:.1f}) + 20*log10({freq_ghz:.1f}) + 92.45 = {fspl_db:.2f} dB\n")

    # --- Step 3: Calculate Received Signal Power at LNA Input ---
    print("Step 3: Calculate Received Signal Power (S) at LNA Input")
    rx_pre_lna_loss_db = rx_ant_loss_db + rx_filter_loss_db
    signal_at_lna_input_dbm = eirp_dbm - fspl_db + rx_ant_gain_db - rx_pre_lna_loss_db
    print(f"S = EIRP (dBm) - FSPL (dB) + Rx Antenna Gain (dB) - Rx Pre-LNA Losses (dB)")
    print(f"S = {eirp_dbm:.2f} - {fspl_db:.2f} + {rx_ant_gain_db:.1f} - ({rx_ant_loss_db:.1f} + {rx_filter_loss_db:.1f}) = {signal_at_lna_input_dbm:.2f} dBm\n")

    # --- Step 4: Calculate Receiver Cascaded Noise Figure (from LNA onwards) ---
    print("Step 4: Calculate Receiver Cascaded Noise Figure (NF)")
    g1_lin = db_to_linear(lna_gain_db)
    f1_lin = db_to_linear(lna_nf_db)
    g2_lin = db_to_linear(-mixer_loss_db)
    f2_lin = db_to_linear(mixer_loss_db) # NF of passive mixer is its loss
    g3_lin = db_to_linear(-if_filter_loss_db)
    f3_lin = db_to_linear(if_filter_loss_db) # NF of passive filter is its loss
    g4_lin = db_to_linear(if_amp_gain_db)
    f4_lin = db_to_linear(if_amp_nf_db)

    # Friis' formula for cascaded noise figure: F = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
    nf_cascade_lin = f1_lin + (f2_lin - 1)/g1_lin + (f3_lin - 1)/(g1_lin * g2_lin) + (f4_lin - 1)/(g1_lin * g2_lin * g3_lin)
    nf_cascade_db = linear_to_db(nf_cascade_lin)
    print("Using Friis' formula for the chain starting from the LNA.")
    print(f"Receiver Noise Figure (from LNA onwards) = {nf_cascade_db:.2f} dB\n")

    # --- Step 5: Calculate Thermal Noise Power ---
    print("Step 5: Calculate Thermal Noise Power (N_thermal)")
    bw_hz = bw_khz * 1000
    thermal_noise_dbm = 10 * math.log10(k_boltzmann * temp_k * bw_hz * 1000)
    print(f"N_thermal = 10*log10(k * T * BW * 1000)")
    print(f"N_thermal = 10*log10({k_boltzmann:.2e} * {temp_k:.1f} * {bw_hz:.0f} * 1000) = {thermal_noise_dbm:.2f} dBm\n")

    # --- Step 6: Calculate Total Noise Power at LNA Input ---
    print("Step 6: Calculate Total Noise Power (N) at LNA Input")
    total_noise_at_lna_input_dbm = thermal_noise_dbm + nf_cascade_db
    print(f"N = Thermal Noise (dBm) + Receiver NF (dB)")
    print(f"N = {thermal_noise_dbm:.2f} + {nf_cascade_db:.2f} = {total_noise_at_lna_input_dbm:.2f} dBm\n")

    # --- Step 7: Calculate Final SNR ---
    print("--- FINAL RESULT ---")
    final_snr_db = signal_at_lna_input_dbm - total_noise_at_lna_input_dbm
    print(f"Final SNR = Signal Power (dBm) - Total Noise Power (dBm)")
    print(f"Final SNR = {signal_at_lna_input_dbm:.2f} - ({total_noise_at_lna_input_dbm:.2f}) = {final_snr_db:.2f} dB")
    
    return final_snr_db

if __name__ == '__main__':
    snr = solve_snr()
    # The final answer in the required format.
    # print(f"\n<<<{snr:.2f}>>>")