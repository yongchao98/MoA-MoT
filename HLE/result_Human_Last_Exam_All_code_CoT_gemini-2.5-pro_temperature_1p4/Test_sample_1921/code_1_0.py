import math

def solve_snr():
    """
    Calculates the Signal-to-Noise Ratio (SNR) for the given microwave link budget.
    """

    # --- Given Parameters ---
    # Tx
    p_tx_dbm = 30.0
    g_tx_ant_db = 20.0
    l_tx_ant_db = 1.0
    l_tx_filter_db = 1.0
    l_tx_cable_db = 1.0
    freq_hz = 24e9
    bw_hz = 100e3

    # Path
    distance_km = 10.0

    # Rx
    g_rx_ant_db = 1.0
    l_rx_ant_db = 0.5
    l_rx_filter_db = 1.0
    lna_gain_db = 36.0
    lna_nf_db = 2.0
    mixer_loss_db = 9.0
    # Mixer NF is assumed equal to its conversion loss
    mixer_nf_db = 9.0
    if_filter_loss_db = 1.0
    # Passive filter NF is equal to its loss
    if_filter_nf_db = 1.0
    if_amp_gain_db = 23.0
    if_amp_nf_db = 0.0  # Negligible NF

    # Environment
    temp_k = 300.0

    # --- Constants ---
    k_boltzmann = 1.38e-23

    # --- Helper Functions ---
    def to_linear(db_val): return 10**(db_val / 10)
    def to_db(lin_val): return 10 * math.log10(lin_val) if lin_val > 0 else -float('inf')


    # --- Step-by-step Calculation ---
    print("Calculating Signal-to-Noise Ratio (SNR)\n")

    # 1. Effective Isotropic Radiated Power (EIRP)
    l_tx_total_db = l_tx_ant_db + l_tx_filter_db + l_tx_cable_db
    eirp_dbm = p_tx_dbm + g_tx_ant_db - l_tx_total_db
    print("Step 1: Calculate Transmitter EIRP")
    print("EIRP (dBm) = Tx Power (dBm) + Tx Ant Gain (dB) - Tx Losses (dB)")
    print(f"EIRP (dBm) = {p_tx_dbm} + {g_tx_ant_db} - ({l_tx_ant_db} + {l_tx_filter_db} + {l_tx_cable_db}) = {eirp_dbm:.2f} dBm\n")

    # 2. Free Space Path Loss (FSPL)
    freq_mhz = freq_hz / 1e6
    fspl_db = 32.45 + 20 * math.log10(freq_mhz) + 20 * math.log10(distance_km)
    print("Step 2: Calculate Free Space Path Loss (FSPL)")
    print("FSPL (dB) = 32.45 + 20*log10(Freq_MHz) + 20*log10(Dist_km)")
    print(f"FSPL (dB) = 32.45 + 20*log10({freq_mhz:.0f}) + 20*log10({distance_km}) = {fspl_db:.2f} dB\n")

    # 3. Received Signal Power at Rx Antenna
    p_rx_dbm = eirp_dbm - fspl_db + g_rx_ant_db
    print("Step 3: Calculate Received Power at Rx Antenna")
    print("Rx Power (dBm) = EIRP (dBm) - FSPL (dB) + Rx Ant Gain (dB)")
    print(f"Rx Power (dBm) = {eirp_dbm:.2f} - {fspl_db:.2f} + {g_rx_ant_db} = {p_rx_dbm:.2f} dBm\n")

    # 4. Input-Referred Noise Power
    pn_watts = k_boltzmann * temp_k * bw_hz
    pn_dbm = to_db(pn_watts * 1000) # convert W to mW for dBm
    print("Step 4: Calculate Input Thermal Noise Power (kTB)")
    print("Noise Power (dBm) = 10*log10(k * T * B * 1000)")
    print(f"Noise Power (dBm) = 10*log10({k_boltzmann:.2e} * {temp_k} * {bw_hz} * 1000) = {pn_dbm:.2f} dBm\n")

    # 5. Input SNR
    snr_in_db = p_rx_dbm - pn_dbm
    print("Step 5: Calculate Input SNR")
    print("Input SNR (dB) = Rx Power (dBm) - Noise Power (dBm)")
    print(f"Input SNR (dB) = {p_rx_dbm:.2f} - ({pn_dbm:.2f}) = {snr_in_db:.2f} dB\n")

    # 6. Total System Noise Figure (NF)
    lna_gain_lin = to_linear(lna_gain_db)
    lna_nf_lin = to_linear(lna_nf_db)
    mixer_gain_lin = to_linear(-mixer_loss_db)
    mixer_nf_lin = to_linear(mixer_nf_db)
    if_filter_gain_lin = to_linear(-if_filter_loss_db)
    if_filter_nf_lin = to_linear(if_filter_nf_db)

    # Friis formula for the cascade starting from the LNA
    f_cascade_lin = lna_nf_lin + (mixer_nf_lin - 1) / lna_gain_lin + \
                    (if_filter_nf_lin - 1) / (lna_gain_lin * mixer_gain_lin)
    nf_cascade_db = to_db(f_cascade_lin)

    # Total NF is the sum of pre-LNA losses and the cascaded NF
    l_pre_lna_db = l_rx_ant_db + l_rx_filter_db
    nf_total_db = l_pre_lna_db + nf_cascade_db
    print("Step 6: Calculate Total System Noise Figure (NF)")
    print("NF is sum of losses before LNA plus the cascaded NF of active components.")
    print(f"Pre-LNA Losses (dB) = Ant Loss + Filter Loss = {l_rx_ant_db} + {l_rx_filter_db} = {l_pre_lna_db:.2f} dB")
    print(f"LNA+ Chain NF (dB) = {nf_cascade_db:.2f} dB (from Friis formula)")
    print(f"Total System NF (dB) = {l_pre_lna_db:.2f} + {nf_cascade_db:.2f} = {nf_total_db:.2f} dB\n")

    # 7. Final Output SNR
    snr_out_db = snr_in_db - nf_total_db
    print("Step 7: Calculate Final Output SNR")
    print("Output SNR (dB) = Input SNR (dB) - Total System NF (dB)")
    print(f"Output SNR (dB) = {snr_in_db:.2f} - {nf_total_db:.2f} = {snr_out_db:.2f} dB")
    
    return snr_out_db

if __name__ == '__main__':
    final_snr = solve_snr()
    # The final answer is wrapped in <<<>>>
    # print(f"\n<<<{final_snr:.2f}>>>")