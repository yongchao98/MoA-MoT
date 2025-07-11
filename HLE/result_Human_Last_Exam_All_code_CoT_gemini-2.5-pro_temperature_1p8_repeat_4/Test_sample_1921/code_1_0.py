import math

def run_calculation():
    """
    Calculates the SNR of a microwave link based on provided parameters.
    """

    # --- Given Parameters ---
    # Transmitter (Tx)
    tx_power_dbm = 30.0  # dBm
    tx_ant_gain_db = 20.0  # dB
    tx_ant_loss_db = 1.0  # dB
    tx_filter_loss_db = 1.0  # dB
    tx_cable_loss_db = 1.0  # dB
    
    # Path
    freq_ghz = 24.0  # GHz
    dist_km = 10.0  # km

    # Receiver (Rx)
    rx_ant_gain_db = 1.0  # dB
    rx_ant_loss_db = 0.5  # dB
    rx_filter1_loss_db = 1.0 # dB
    lna_gain_db = 36.0  # dB
    lna_nf_db = 2.0  # dB
    mixer_loss_db = 9.0  # dB
    if_filter_loss_db = 1.0 # dB
    if_amp_gain_db = 23.0 # dB
    if_amp_nf_db = 0.0 # dB (negligible)

    # System
    temp_k = 300.0  # K
    bandwidth_hz = 100_000.0  # 100 kHz

    # --- Constants ---
    BOLTZMANN_K = 1.380649e-23  # J/K

    print("--- Step 1: Calculate Effective Isotropic Radiated Power (EIRP) ---")
    tx_total_loss_db = tx_ant_loss_db + tx_filter_loss_db + tx_cable_loss_db
    eirp_dbm = tx_power_dbm + tx_ant_gain_db - tx_total_loss_db
    print(f"Tx Losses = {tx_ant_loss_db} dB + {tx_filter_loss_db} dB + {tx_cable_loss_db} dB = {tx_total_loss_db:.1f} dB")
    print(f"EIRP (dBm) = Tx Power + Tx Antenna Gain - Tx Losses")
    print(f"EIRP (dBm) = {tx_power_dbm} dBm + {tx_ant_gain_db} dB - {tx_total_loss_db:.1f} dB = {eirp_dbm:.2f} dBm\n")

    print("--- Step 2: Calculate Free Space Path Loss (FSPL) ---")
    fspl_db = 20 * math.log10(dist_km) + 20 * math.log10(freq_ghz) + 92.45
    print(f"FSPL (dB) = 20*log10(distance_km) + 20*log10(frequency_GHz) + 92.45")
    print(f"FSPL (dB) = 20*log10({dist_km}) + 20*log10({freq_ghz}) + 92.45 = {fspl_db:.2f} dB\n")

    print("--- Step 3: Calculate Received Signal Power (S) ---")
    prx_dbm = eirp_dbm - fspl_db + rx_ant_gain_db
    print("Received Power S (dBm) = EIRP - FSPL + Rx Antenna Gain")
    print(f"S (dBm) = {eirp_dbm:.2f} dBm - {fspl_db:.2f} dB + {rx_ant_gain_db} dB = {prx_dbm:.2f} dBm\n")

    print("--- Step 4: Calculate Total Receiver Noise Figure (NF) ---")
    # Helper functions for conversion
    nf_db_to_linear = lambda nf: 10**(nf / 10)
    gain_db_to_linear = lambda g: 10**(g / 10)

    # List of components in the receiver chain for NF calculation
    # Each component is a tuple: (Name, Gain_dB, NF_dB)
    # Note: Loss is negative gain, and NF of a passive component is its loss.
    rx_chain = [
        ("Rx Ant Loss", -rx_ant_loss_db, rx_ant_loss_db),
        ("Rx Input Filter", -rx_filter1_loss_db, rx_filter1_loss_db),
        ("LNA", lna_gain_db, lna_nf_db),
        ("Mixer", -mixer_loss_db, mixer_loss_db),
        ("IF Filter", -if_filter_loss_db, if_filter_loss_db),
        ("IF Amplifier", if_amp_gain_db, if_amp_nf_db),
    ]

    # Friis formula calculation
    cascaded_f = 0
    cascaded_g = 1
    print("Calculating cascaded noise figure using Friis formula F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...")
    for i, (name, g_db, nf_db) in enumerate(rx_chain):
        f_lin = nf_db_to_linear(nf_db)
        term = (f_lin - 1) / cascaded_g if i > 0 else f_lin
        cascaded_f += term
        
        g_lin = gain_db_to_linear(g_db)
        cascaded_g *= g_lin
        
    nf_total_db = 10 * math.log10(cascaded_f)
    print(f"Total Receiver Noise Figure = {nf_total_db:.2f} dB\n")

    print("--- Step 5: Calculate Total Input-Referred Noise Power (N) ---")
    noise_thermal_watts = BOLTZMANN_K * temp_k * bandwidth_hz
    noise_thermal_dbm = 10 * math.log10(noise_thermal_watts * 1000)
    print(f"Thermal Noise (dBm) = 10*log10(k * T * B * 1000)")
    print(f"Thermal Noise (dBm) = 10*log10({BOLTZMANN_K:.2e} * {temp_k} * {bandwidth_hz} * 1000) = {noise_thermal_dbm:.2f} dBm")
    
    total_noise_dbm = noise_thermal_dbm + nf_total_db
    print("Total Noise N (dBm) = Thermal Noise + Total Receiver NF")
    print(f"N (dBm) = {noise_thermal_dbm:.2f} dBm + {nf_total_db:.2f} dB = {total_noise_dbm:.2f} dBm\n")

    print("--- Step 6: Calculate Final Signal-to-Noise Ratio (SNR) ---")
    snr_db = prx_dbm - total_noise_dbm
    print("SNR (dB) = Signal Power S (dBm) - Total Noise N (dBm)")
    print(f"SNR (dB) = {prx_dbm:.2f} dBm - ({total_noise_dbm:.2f} dBm) = {snr_db:.2f} dB")
    
    return snr_db

if __name__ == "__main__":
    final_snr = run_calculation()
    # The final answer is wrapped in <<<>>> as requested.
    print(f"\n<<<SNR = {final_snr:.2f} dB>>>")
