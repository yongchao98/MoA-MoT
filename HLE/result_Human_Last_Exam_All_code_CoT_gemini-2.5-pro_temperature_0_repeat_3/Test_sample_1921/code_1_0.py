import math

def main():
    """
    Calculates the Signal-to-Noise Ratio (SNR) for a microwave link based on given parameters.
    """
    # --- GIVEN PARAMETERS ---

    # Constants
    k = 1.380649e-23  # Boltzmann constant (J/K)
    T = 300           # Ambient temperature (K)
    c = 299792458     # Speed of light (m/s)

    # Transmitter (Tx)
    tx_power_dbm = 30.0
    tx_ant_gain_db = 20.0
    tx_ant_loss_db = 1.0
    tx_filter_loss_db = 1.0
    tx_cable_loss_db = 1.0

    # Path
    freq_hz = 24e9  # 24 GHz
    distance_m = 10e3 # 10 km
    bandwidth_hz = 100e3 # 100 kHz

    # Receiver (Rx)
    rx_ant_gain_db = 1.0
    rx_ant_loss_db = 0.5
    rx_filter_loss_db = 1.0
    lna_gain_db = 36.0
    lna_nf_db = 2.0
    mixer_loss_db = 9.0
    mixer_nf_db = 9.0  # NF of a passive mixer equals its conversion loss
    if_filter_loss_db = 1.0
    if_amp_gain_db = 23.0
    if_amp_nf_db = 0.0  # Negligible NF
    output_filter_loss_db = 1.0

    # --- CALCULATIONS ---

    # 1. Calculate EIRP (Effective Isotropic Radiated Power)
    total_tx_loss_db = tx_ant_loss_db + tx_filter_loss_db + tx_cable_loss_db
    eirp_dbm = tx_power_dbm + tx_ant_gain_db - total_tx_loss_db

    # 2. Calculate FSPL (Free Space Path Loss)
    # FSPL (dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
    fspl_db = 20 * math.log10(distance_m) + 20 * math.log10(freq_hz) + 20 * math.log10(4 * math.pi / c)

    # 3. Calculate Received Signal Power (at antenna input, before Rx losses)
    prx_dbm = eirp_dbm - fspl_db + rx_ant_gain_db

    # 4. Calculate Receiver System Noise Figure (NF_sys) using Friis formula
    # F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
    # where F is noise factor (linear) and G is gain (linear)
    
    # List of components in the receiver chain: (Gain in dB, Noise Figure in dB)
    # A lossy component has a gain < 0 dB and a noise figure equal to its loss.
    rx_chain = [
        (-rx_ant_loss_db, rx_ant_loss_db),
        (-rx_filter_loss_db, rx_filter_loss_db),
        (lna_gain_db, lna_nf_db),
        (-mixer_loss_db, mixer_nf_db),
        (-if_filter_loss_db, if_filter_loss_db),
        (if_amp_gain_db, if_amp_nf_db),
        (-output_filter_loss_db, output_filter_loss_db)
    ]

    total_f = 0
    cumulative_g = 1.0
    for g_db, nf_db in rx_chain:
        g_lin = 10**(g_db / 10)
        f_lin = 10**(nf_db / 10)
        
        if total_f == 0: # First stage
            total_f = f_lin
        else:
            total_f += (f_lin - 1) / cumulative_g
        
        cumulative_g *= g_lin
        
    nf_sys_db = 10 * math.log10(total_f)

    # 5. Calculate Thermal Noise Power
    # N = k * T * B
    thermal_noise_dbm = 10 * math.log10(k * T * bandwidth_hz * 1000) # *1000 to convert W to mW

    # 6. Calculate Total Input-Referred Noise Power
    n_total_dbm = thermal_noise_dbm + nf_sys_db

    # 7. Calculate Final SNR
    snr_db = prx_dbm - n_total_dbm

    # --- OUTPUT RESULTS ---
    print("--- Link Budget Analysis ---")
    print(f"\n1. Effective Isotropic Radiated Power (EIRP)")
    print(f"EIRP = Tx Power + Tx Ant Gain - Tx Losses")
    print(f"EIRP = {tx_power_dbm:.1f} dBm + {tx_ant_gain_db:.1f} dB - {total_tx_loss_db:.1f} dB = {eirp_dbm:.2f} dBm")

    print(f"\n2. Free Space Path Loss (FSPL)")
    print(f"FSPL for {distance_m/1000:.1f} km at {freq_hz/1e9:.1f} GHz = {fspl_db:.2f} dB")

    print(f"\n3. Received Signal Power (Prx)")
    print(f"Prx = EIRP - FSPL + Rx Ant Gain")
    print(f"Prx = {eirp_dbm:.2f} dBm - {fspl_db:.2f} dB + {rx_ant_gain_db:.1f} dB = {prx_dbm:.2f} dBm")

    print(f"\n4. Total Input-Referred Noise Power (N_total)")
    print(f"Thermal Noise in {bandwidth_hz/1000:.1f} kHz BW = {thermal_noise_dbm:.2f} dBm")
    print(f"Receiver System Noise Figure = {nf_sys_db:.2f} dB")
    print(f"N_total = Thermal Noise + System NF")
    print(f"N_total = {thermal_noise_dbm:.2f} dBm + {nf_sys_db:.2f} dB = {n_total_dbm:.2f} dBm")

    print(f"\n--- Final SNR Calculation ---")
    print(f"SNR = Received Signal Power - Total Input-Referred Noise Power")
    print(f"SNR = {prx_dbm:.2f} dBm - ({n_total_dbm:.2f} dBm)")
    print(f"Resulting SNR = {snr_db:.2f} dB")

if __name__ == "__main__":
    main()
<<<28.27>>>