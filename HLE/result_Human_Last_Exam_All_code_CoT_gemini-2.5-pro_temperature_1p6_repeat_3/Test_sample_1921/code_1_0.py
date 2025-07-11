import math

def calculate_snr():
    """
    Calculates the Signal-to-Noise Ratio (SNR) for a microwave link based on given parameters.
    """
    # --- Constants ---
    K_BOLTZMANN = 1.380649e-23  # Joules per Kelvin
    C_LIGHT = 299792458        # Meters per second

    # --- Input Parameters ---
    # Transmitter (Tx)
    p_tx_dbm = 30.0
    g_ant_tx_db = 20.0
    l_ant_tx_db = 1.0
    l_filter_tx_db = 1.0
    l_cable_tx_db = 1.0

    # Path
    freq_hz = 24e9
    dist_m = 10e3

    # Receiver (Rx)
    g_ant_rx_db = 1.0
    l_ant_rx_db = 0.5
    l_filter_rx_in_db = 1.0
    nf_lna_db = 2.0
    # Note: Components after the high-gain LNA have negligible impact on the system Noise Figure.

    # System
    bandwidth_hz = 100e3
    t_ambient_k = 300.0
    
    print("--- Step 1: Calculate Effective Isotropic Radiated Power (EIRP) ---")
    l_tx_total_db = l_ant_tx_db + l_filter_tx_db + l_cable_tx_db
    eirp_dbm = p_tx_dbm + g_ant_tx_db - l_tx_total_db
    print(f"EIRP = Tx Power + Tx Antenna Gain - Tx Losses")
    print(f"EIRP = {p_tx_dbm} dBm + {g_ant_tx_db} dB - ({l_ant_tx_db} + {l_filter_tx_db} + {l_cable_tx_db}) dB = {eirp_dbm:.2f} dBm\n")

    print("--- Step 2: Calculate Free Space Path Loss (FSPL) ---")
    # FSPL (dB) = 20*log10(d) + 20*log10(f) + 20*log10(4π/c)
    fspl_db = 20 * math.log10(dist_m) + 20 * math.log10(freq_hz) + 20 * math.log10(4 * math.pi / C_LIGHT)
    print(f"FSPL = 20*log10({dist_m:.0f} m) + 20*log10({freq_hz:.0f} Hz) + 20*log10(4π/c)")
    print(f"FSPL = {fspl_db:.2f} dB\n")

    print("--- Step 3: Calculate Received Signal Power (S) at antenna connector ---")
    s_dbm = eirp_dbm - fspl_db + g_ant_rx_db
    print(f"Signal Power = EIRP - FSPL + Rx Antenna Gain")
    print(f"Signal Power = {eirp_dbm:.2f} dBm - {fspl_db:.2f} dB + {g_ant_rx_db} dB = {s_dbm:.2f} dBm\n")

    print("--- Step 4: Calculate System Noise Figure (NF_sys) ---")
    # The system noise figure is dominated by the losses before the LNA and the LNA's own NF.
    # NF_sys = L_ant_rx + L_filter_rx + NF_lna
    nf_sys_db = l_ant_rx_db + l_filter_rx_in_db + nf_lna_db
    print("System NF = Rx Antenna Loss + Rx Filter Loss + LNA NF")
    print(f"System NF = {l_ant_rx_db} dB + {l_filter_rx_in_db} dB + {nf_lna_db} dB = {nf_sys_db:.2f} dB\n")

    print("--- Step 5: Calculate Thermal Noise Power (N_th) ---")
    n_th_watts = K_BOLTZMANN * t_ambient_k * bandwidth_hz
    n_th_dbm = 10 * math.log10(n_th_watts / 0.001)
    print(f"Thermal Noise = 10 * log10(k * T * B / 1mW)")
    print(f"Thermal Noise = 10 * log10({K_BOLTZMANN:.3e} * {t_ambient_k}K * {bandwidth_hz:.0f}Hz / 0.001) = {n_th_dbm:.2f} dBm\n")

    print("--- Step 6: Calculate Total Input-Referred Noise (N) ---")
    n_total_dbm = n_th_dbm + nf_sys_db
    print("Total Noise = Thermal Noise + System NF")
    print(f"Total Noise = {n_th_dbm:.2f} dBm + {nf_sys_db:.2f} dB = {n_total_dbm:.2f} dBm\n")
    
    print("--- Step 7: Calculate Final Signal-to-Noise Ratio (SNR) ---")
    snr_db = s_dbm - n_total_dbm
    print("SNR = Signal Power - Total Noise")
    print(f"Final SNR = {s_dbm:.2f} dBm - ({n_total_dbm:.2f} dBm) = {snr_db:.2f} dB")
    
    return snr_db

if __name__ == '__main__':
    final_snr = calculate_snr()
    # The final answer is wrapped in <<<>>>
    print(f"\n<<<{final_snr:.2f}>>>")