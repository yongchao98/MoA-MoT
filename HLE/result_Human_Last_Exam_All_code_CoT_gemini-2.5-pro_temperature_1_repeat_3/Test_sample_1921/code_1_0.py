import math

def solve_snr():
    """
    Calculates the SNR of a microwave link based on the provided parameters.
    """
    # --- Constants ---
    k = 1.380649e-23  # Boltzmann's constant (J/K)
    c = 2.998e8       # Speed of light (m/s)

    # --- Given Parameters ---
    # Transmitter (Tx)
    P_tx_dbm = 30.0       # Tx Power (dBm)
    G_tx_db = 20.0        # Tx Antenna Gain (dB)
    L_ant_tx_db = 1.0     # Tx Antenna Loss (dB)
    L_filter_tx_db = 1.0  # Tx Output Filter Loss (dB)
    L_cable_tx_db = 1.0   # Tx Cable Loss (dB)
    freq_hz = 24e9        # Frequency (Hz)
    BW_hz = 100e3         # Signal Bandwidth (Hz)

    # Path
    dist_m = 10e3         # Distance (m)

    # Receiver (Rx)
    G_rx_db = 1.0         # Rx Antenna Gain (dB)
    L_ant_rx_db = 0.5     # Rx Antenna Loss (dB)
    L_filter1_rx_db = 1.0 # Rx Input Filter Loss (dB)
    G_lna_db = 36.0       # LNA Gain (dB)
    NF_lna_db = 2.0       # LNA Noise Figure (dB)
    L_mixer_db = 9.0      # Mixer Conversion Loss (dB)
    NF_mixer_db = 9.0     # Mixer Noise Figure (dB)
    L_filter2_rx_db = 1.0 # IF Filter Loss (dB)
    G_if_db = 23.0        # IF Amplifier Gain (dB)
    NF_if_db = 0.0        # IF Amplifier Noise Figure (dB, negligible)
    # The final output filter does not affect the SNR calculation at the IF stage.

    # Environment
    T_kelvin = 300.0      # Ambient Temperature (K)

    # --- Step 1: Calculate Transmitter EIRP ---
    L_tx_total_db = L_filter_tx_db + L_cable_tx_db + L_ant_tx_db
    EIRP_dbm = P_tx_dbm + G_tx_db - L_tx_total_db

    # --- Step 2: Calculate Free Space Path Loss (FSPL) ---
    # FSPL (dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
    fspl_db = 20 * math.log10(dist_m) + 20 * math.log10(freq_hz) + 20 * math.log10(4 * math.pi / c)

    # --- Step 3: Calculate Received Signal Power (S) ---
    # Reference point is at the input to the receiver chain (after antenna gain, before other Rx losses)
    S_dbm = EIRP_dbm - fspl_db + G_rx_db

    # --- Step 4: Calculate Total Receiver Noise Figure (NF_rx) ---
    # The total NF is calculated at the same reference point as the signal power.
    # We use the Friis formula for cascaded noise figure: F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
    # where F and G are linear values (not dB).

    # Convert dB values to linear for the components in the Rx chain
    # Stage 1: Rx Antenna Loss
    L1_lin = 10**(L_ant_rx_db / 10)
    G1_lin = 1 / L1_lin
    F1_lin = L1_lin # Noise Factor of a passive lossy component is its loss

    # Stage 2: Rx Input Filter Loss
    L2_lin = 10**(L_filter1_rx_db / 10)
    G2_lin = 1 / L2_lin
    F2_lin = L2_lin

    # Stage 3: LNA
    G3_lin = 10**(G_lna_db / 10)
    F3_lin = 10**(NF_lna_db / 10)

    # Stage 4: Mixer
    G4_lin = 1 / (10**(L_mixer_db / 10))
    F4_lin = 10**(NF_mixer_db / 10)
    
    # Later stages will have negligible contribution due to high LNA gain, but we can include them
    # Stage 5: IF Filter
    L5_lin = 10**(L_filter2_rx_db / 10)
    G5_lin = 1 / L5_lin
    F5_lin = L5_lin

    # Stage 6: IF Amplifier
    G6_lin = 10**(G_if_db / 10)
    F6_lin = 10**(NF_if_db / 10)

    # Apply Friis formula
    F_total_lin = F1_lin + (F2_lin - 1)/G1_lin + \
                  (F3_lin - 1)/(G1_lin * G2_lin) + \
                  (F4_lin - 1)/(G1_lin * G2_lin * G3_lin) + \
                  (F5_lin - 1)/(G1_lin * G2_lin * G3_lin * G4_lin) + \
                  (F6_lin - 1)/(G1_lin * G2_lin * G3_lin * G4_lin * G5_lin)

    NF_rx_total_db = 10 * math.log10(F_total_lin)

    # --- Step 5: Calculate Total Noise Power (N) ---
    # Thermal noise in dBm
    N_thermal_watts = k * T_kelvin * BW_hz
    N_thermal_dbm = 10 * math.log10(N_thermal_watts / 0.001)

    # Total noise power at the reference point
    N_total_dbm = N_thermal_dbm + NF_rx_total_db

    # --- Step 6: Calculate Final SNR ---
    SNR_db = S_dbm - N_total_dbm

    # --- Print the results step-by-step ---
    print("--- Link Budget Calculation ---")
    print("\n1. Transmitter and Path Loss Calculation:")
    print(f"Tx Power = {P_tx_dbm:.2f} dBm")
    print(f"Tx Antenna Gain = {G_tx_db:.2f} dB")
    print(f"Tx Total Loss = {L_tx_total_db:.2f} dB")
    print(f"Effective Isotropic Radiated Power (EIRP) = {P_tx_dbm:.2f} + {G_tx_db:.2f} - {L_tx_total_db:.2f} = {EIRP_dbm:.2f} dBm")
    print(f"Free Space Path Loss (FSPL) = {fspl_db:.2f} dB")

    print("\n2. Signal Power Calculation:")
    print(f"Rx Antenna Gain = {G_rx_db:.2f} dB")
    print(f"Signal Power (S) = EIRP - FSPL + G_rx")
    print(f"Signal Power (S) = {EIRP_dbm:.2f} dBm - {fspl_db:.2f} dB + {G_rx_db:.2f} dB = {S_dbm:.2f} dBm")

    print("\n3. Noise Power Calculation:")
    print(f"Thermal Noise Power (at {T_kelvin}K, {BW_hz/1000}kHz BW) = {N_thermal_dbm:.2f} dBm")
    print(f"Total Receiver Noise Figure = {NF_rx_total_db:.2f} dB")
    print(f"Total Noise Power (N) = Thermal Noise + Rx Noise Figure")
    print(f"Total Noise Power (N) = {N_thermal_dbm:.2f} dBm + {NF_rx_total_db:.2f} dB = {N_total_dbm:.2f} dBm")

    print("\n4. Final SNR Calculation:")
    print(f"SNR (dB) = Signal Power (dBm) - Noise Power (dBm)")
    print(f"SNR (dB) = {S_dbm:.2f} - ({N_total_dbm:.2f})")
    print(f"Final SNR = {SNR_db:.2f} dB")
    
    # Final answer in the required format
    print(f"\n<<<{SNR_db:.2f}>>>")

solve_snr()