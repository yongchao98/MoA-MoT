import math

def calculate_snr():
    """
    Calculates the SNR of a microwave link based on given parameters.
    The function prints each step of the calculation.
    """
    # --- GIVEN PARAMETERS ---
    # Transmitter
    P_tx_dbm = 30.0    # Transmitter Power (dBm)
    G_tx_ant_db = 20.0   # Tx Antenna Gain (dB)
    L_tx_ant_db = 1.0    # Tx Antenna Loss (dB)
    L_tx_filter_db = 1.0 # Tx Filter Loss (dB)
    L_tx_cable_db = 1.0  # Tx Cable Loss (dB)

    # Path
    f_hz = 24e9        # Frequency (Hz)
    d_m = 10e3         # Distance (m)

    # Receiver
    G_rx_ant_db = 1.0    # Rx Antenna Gain (dB)
    L_rx_ant_db = 0.5    # Rx Antenna Loss (dB)
    L_rx_filter_db = 1.0 # Rx Input Filter Loss (dB)
    G_lna_db = 36.0      # LNA Gain (dB)
    NF_lna_db = 2.0      # LNA Noise Figure (dB)
    L_mixer_db = 9.0     # Mixer Conversion Loss (dB)
    NF_mixer_db = 9.0    # Mixer Noise Figure (assumed equal to loss)
    L_if_filter_db = 1.0 # IF Filter Loss (dB)
    NF_if_filter_db = 1.0# IF Filter Noise Figure (assumed equal to loss)
    G_if_amp_db = 23.0   # IF Amplifier Gain (dB)
    NF_if_amp_db = 0.0   # IF Amplifier Noise Figure (negligible)
    
    # System
    T_kelvin = 300     # Ambient Temperature (K)
    B_hz = 100e3       # Signal Bandwidth (Hz)
    
    # --- CONSTANTS ---
    k = 1.380649e-23   # Boltzmann's constant (J/K)
    c = 299792458      # Speed of light (m/s)

    # --- STEP 1: Calculate Effective Isotropic Radiated Power (EIRP) ---
    print("--- 1. Effective Isotropic Radiated Power (EIRP) Calculation ---")
    L_tx_total_db = L_tx_ant_db + L_tx_filter_db + L_tx_cable_db
    eirp_dbm = P_tx_dbm + G_tx_ant_db - L_tx_total_db
    print(f"EIRP (dBm) = Tx Power (dBm) + Tx Antenna Gain (dB) - Total Tx Loss (dB)")
    print(f"EIRP = {P_tx_dbm} + {G_tx_ant_db} - ({L_tx_ant_db} + {L_tx_filter_db} + {L_tx_cable_db})")
    print(f"EIRP = {eirp_dbm:.2f} dBm\n")

    # --- STEP 2: Calculate Free Space Path Loss (FSPL) ---
    print("--- 2. Free Space Path Loss (FSPL) Calculation ---")
    # Using formula: FSPL(dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
    # The constant 20*log10(4*pi/c) is approximately -147.55 dB
    fspl_db = 20 * math.log10(d_m) + 20 * math.log10(f_hz) - 147.55
    print(f"FSPL (dB) = 20*log10(distance_m) + 20*log10(frequency_hz) - 147.55")
    print(f"FSPL = 20*log10({d_m:.0f}) + 20*log10({f_hz:.0f}) - 147.55")
    print(f"FSPL = {fspl_db:.2f} dB\n")

    # --- STEP 3: Calculate Received Signal Power (S) ---
    print("--- 3. Received Signal Power (S) Calculation ---")
    S_dbm = eirp_dbm - fspl_db + G_rx_ant_db
    print("S (dBm) = EIRP (dBm) - FSPL (dB) + Rx Antenna Gain (dB)")
    print(f"S = {eirp_dbm:.2f} - {fspl_db:.2f} + {G_rx_ant_db}")
    print(f"S = {S_dbm:.2f} dBm\n")

    # --- STEP 4: Calculate System Noise Figure (NF_sys) ---
    print("--- 4. System Noise Figure (NF_sys) Calculation ---")
    
    # Convert dB to linear values for Friis formula
    # F = 10^(NF/10), G = 10^(G/10)
    # Note: For passive lossy components, NF = Loss
    nf_rx_ant_lin = 10**(L_rx_ant_db / 10)
    g_rx_ant_lin = 10**(-L_rx_ant_db / 10)
    nf_rx_filter_lin = 10**(L_rx_filter_db / 10)
    g_rx_filter_lin = 10**(-L_rx_filter_db / 10)
    nf_lna_lin = 10**(NF_lna_db / 10)
    g_lna_lin = 10**(G_lna_db / 10)
    nf_mixer_lin = 10**(NF_mixer_db / 10)
    g_mixer_lin = 10**(-L_mixer_db / 10) # Mixer has conversion loss, so negative gain
    nf_if_filter_lin = 10**(NF_if_filter_db / 10)

    # Friis Formula for Noise Factor: F_sys = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
    F_sys_lin = nf_rx_ant_lin + \
                (nf_rx_filter_lin - 1) / g_rx_ant_lin + \
                (nf_lna_lin - 1) / (g_rx_ant_lin * g_rx_filter_lin) + \
                (nf_mixer_lin - 1) / (g_rx_ant_lin * g_rx_filter_lin * g_lna_lin)
                # Subsequent stages have negligible contribution due to high LNA gain
    
    nf_sys_db = 10 * math.log10(F_sys_lin)
    print("Using the Friis formula for cascaded noise figure:")
    print("NF_sys = 10*log10(F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...)")
    print(f"The receiver chain starts with components (Loss/NF): Rx Ant Loss={L_rx_ant_db}dB, Rx Filter={L_rx_filter_db}dB, LNA NF={NF_lna_db}dB, ...")
    print(f"Resulting System Noise Factor (linear) F_sys = {F_sys_lin:.3f}")
    print(f"NF_sys (dB) = 10*log10({F_sys_lin:.3f}) = {nf_sys_db:.2f} dB\n")

    # --- STEP 5: Calculate Total Input Noise Power (N) ---
    print("--- 5. Total Input Noise Power (N) Calculation ---")
    N_floor_watts = k * T_kelvin * B_hz
    N_floor_dbm = 10 * math.log10(N_floor_watts / 0.001)
    
    N_total_dbm = N_floor_dbm + nf_sys_db
    
    print("N (dBm) = Thermal Noise (dBm) + System Noise Figure (dB)")
    print(f"Thermal Noise (dBm) = 10*log10(k * T * B * 1000) = 10*log10({k:.2e} * {T_kelvin} * {B_hz:.0f} * 1000) = {N_floor_dbm:.2f} dBm")
    print(f"N = {N_floor_dbm:.2f} + {nf_sys_db:.2f}")
    print(f"N = {N_total_dbm:.2f} dBm\n")

    # --- STEP 6: Calculate Final SNR ---
    print("--- 6. Final Signal-to-Noise Ratio (SNR) Calculation ---")
    snr_db = S_dbm - N_total_dbm
    print("SNR (dB) = S (dBm) - N (dBm)")
    print(f"SNR = {S_dbm:.2f} - ({N_total_dbm:.2f})")
    print(f"SNR = {snr_db:.2f} dB")
    
    return snr_db

if __name__ == '__main__':
    final_snr = calculate_snr()
    # The final answer is wrapped in <<<>>> as requested.
    # print(f"\n<<<Result>>>\n{final_snr:.2f}")
    # Final answer requested in a very specific format
    # The line below will be captured as the final answer.
