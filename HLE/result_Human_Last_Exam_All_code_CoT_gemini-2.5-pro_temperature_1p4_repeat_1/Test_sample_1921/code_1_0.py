import math

def solve_snr():
    """
    Calculates the SNR of a microwave link based on provided parameters.
    """

    # --- Constants ---
    k = 1.380649e-23  # Boltzmann's constant (J/K)
    c = 299792458     # Speed of light (m/s)
    T = 300           # Ambient temperature (K)

    # --- Given Parameters ---
    # Transmitter (Tx)
    tx_power_dbm = 30.0
    tx_ant_gain_db = 20.0
    tx_ant_loss_db = 1.0
    tx_filter_loss_db = 1.0
    tx_cable_loss_db = 1.0
    freq_hz = 24e9
    bandwidth_hz = 100e3

    # Path
    dist_m = 10e3

    # Receiver (Rx)
    rx_ant_gain_db = 1.0
    
    # Receiver chain components for Noise Figure calculation
    # Format: (Name, Gain_dB, Noise_Figure_dB)
    # Note: For passive components, gain is negative loss and NF is equal to loss.
    rx_chain = [
        ('Rx Antenna Loss', -0.5, 0.5),
        ('Rx Input Filter', -1.0, 1.0),
        ('LNA', 36.0, 2.0),
        ('Mixer', -9.0, 9.0),
        ('IF Filter', -1.0, 1.0),
        ('IF Amplifier', 23.0, 0.0), # Negligible NF
        ('Output Filter', -1.0, 1.0),
    ]

    # --- Helper Functions ---
    def db_to_linear(db):
        return 10**(db / 10)

    def linear_to_db(linear):
        return 10 * math.log10(linear)

    print("--- Step 1: Calculate Effective Isotropically Radiated Power (EIRP) ---")
    tx_total_loss_db = tx_ant_loss_db + tx_filter_loss_db + tx_cable_loss_db
    eirp_dbm = tx_power_dbm + tx_ant_gain_db - tx_total_loss_db
    print(f"EIRP (dBm) = Tx Power (dBm) + Tx Ant Gain (dB) - Tx Losses (dB)")
    print(f"EIRP (dBm) = {tx_power_dbm} + {tx_ant_gain_db} - ({tx_ant_loss_db} + {tx_filter_loss_db} + {tx_cable_loss_db}) = {eirp_dbm:.2f} dBm\n")

    print("--- Step 2: Calculate Free Space Path Loss (FSPL) ---")
    # FSPL formula: 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
    fspl_db = 20 * math.log10(dist_m) + 20 * math.log10(freq_hz) + 20 * math.log10(4 * math.pi / c)
    print(f"FSPL (dB) = 20*log10({dist_m:.0f} m) + 20*log10({freq_hz:.0f} Hz) + 20*log10(4*pi/c)")
    print(f"FSPL (dB) = {20 * math.log10(dist_m):.2f} + {20 * math.log10(freq_hz):.2f} - 147.55 = {fspl_db:.2f} dB\n")

    print("--- Step 3: Calculate Received Signal Power (S) at Antenna Terminal ---")
    rx_signal_power_dbm = eirp_dbm - fspl_db + rx_ant_gain_db
    print(f"S (dBm) = EIRP (dBm) - FSPL (dB) + Rx Ant Gain (dB)")
    print(f"S (dBm) = {eirp_dbm:.2f} - {fspl_db:.2f} + {rx_ant_gain_db} = {rx_signal_power_dbm:.2f} dBm\n")

    print("--- Step 4: Calculate Receiver System Noise Figure (NF_sys) ---")
    # Using Friis formula for cascaded noise figure: F_sys = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
    # Reference plane is the input to the first component in the chain (Rx Antenna Loss element)
    # which corresponds to the antenna terminal.
    
    f_sys_linear = 0.0
    g_cascade_linear = 1.0
    
    # Handle the first stage separately
    g1_db, nf1_db = rx_chain[0][1], rx_chain[0][2]
    f1_linear = db_to_linear(nf1_db)
    g1_linear = db_to_linear(g1_db)
    
    f_sys_linear = f1_linear
    g_cascade_linear = g1_linear

    # Loop through the rest of the stages
    for name, g_db, nf_db in rx_chain[1:]:
        f_linear = db_to_linear(nf_db)
        g_linear = db_to_linear(g_db)
        
        f_sys_linear += (f_linear - 1) / g_cascade_linear
        g_cascade_linear *= g_linear

    nf_sys_db = linear_to_db(f_sys_linear)
    print(f"Using Friis formula for the receiver chain.")
    print(f"System Noise Factor F_sys = {f_sys_linear:.3f}")
    print(f"System Noise Figure NF_sys (dB) = 10*log10({f_sys_linear:.3f}) = {nf_sys_db:.2f} dB\n")

    print("--- Step 5: Calculate Total Noise Power (N) ---")
    # N = k * T * B * F_sys
    noise_power_watts = k * T * bandwidth_hz * f_sys_linear
    noise_power_dbm = 10 * math.log10(noise_power_watts / 0.001)
    print(f"N (dBm) = 10*log10(k * T * B * F_sys / 1mW)")
    print(f"N (dBm) = 10*log10({k:.3e} * {T} * {bandwidth_hz:.0f} * {f_sys_linear:.3f} / 0.001) = {noise_power_dbm:.2f} dBm\n")

    print("--- Step 6: Calculate Final Signal-to-Noise Ratio (SNR) ---")
    snr_db = rx_signal_power_dbm - noise_power_dbm
    print(f"SNR (dB) = S (dBm) - N (dBm)")
    print(f"SNR (dB) = {rx_signal_power_dbm:.2f} - ({noise_power_dbm:.2f}) = {snr_db:.2f} dB")
    
    return snr_db

if __name__ == '__main__':
    final_snr = solve_snr()
    print(f"\n<<<Result>>>")
    print(f"The final SNR is {final_snr:.2f} dB.")
    print(f"<<<{final_snr:.2f}>>>")
