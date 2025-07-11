import math

def db_to_linear(db):
    """Converts a value from dB to linear scale."""
    return 10**(db / 10)

def linear_to_db(linear):
    """Converts a value from linear scale to dB."""
    return 10 * math.log10(linear)

def calculate_snr():
    """
    Calculates the SNR of the microwave link based on the provided parameters.
    """
    # --- Given Parameters ---
    # Constants
    k = 1.380649e-23  # Boltzmann constant in J/K
    T = 300           # Ambient temperature in Kelvin
    c = 299792458     # Speed of light in m/s

    # Transmitter (Tx)
    P_tx_dBm = 30.0
    G_ant_tx_dB = 20.0
    L_ant_tx_dB = 1.0
    L_filter_tx_dB = 1.0
    L_cable_tx_dB = 1.0
    f_Hz = 24e9       # Frequency in Hz
    B_Hz = 100e3      # Signal bandwidth in Hz

    # Path
    d_m = 10e3        # Distance in meters

    # Receiver (Rx)
    G_ant_rx_dB = 1.0
    L_ant_rx_dB = 0.5
    L_filter_in_dB = 1.0
    LNA_gain_dB = 36.0
    LNA_NF_dB = 2.0
    Mixer_loss_dB = 9.0
    IF_filter_loss_dB = 1.0
    IF_amp_gain_dB = 23.0
    IF_amp_NF_dB = 0.0 # Negligible NF
    Output_filter_loss_dB = 1.0

    print("--- Step 1: Calculate Effective Isotropic Radiated Power (EIRP) ---")
    # Power at the antenna input after losses
    power_before_ant = P_tx_dBm - L_cable_tx_dB - L_filter_tx_dB
    # EIRP includes antenna gain but subtracts antenna's own physical loss
    eirp_dBm = power_before_ant + G_ant_tx_dB - L_ant_tx_dB
    print(f"EIRP (dBm) = Tx Power - Cable Loss - Filter Loss + Antenna Gain - Antenna Loss")
    print(f"EIRP (dBm) = {P_tx_dBm} dBm - {L_cable_tx_dB} dB - {L_filter_tx_dB} dB + {G_ant_tx_dB} dBi - {L_ant_tx_dB} dB = {eirp_dBm:.2f} dBm\n")

    print("--- Step 2: Calculate Free Space Path Loss (FSPL) ---")
    fspl_dB = 20 * math.log10(d_m) + 20 * math.log10(f_Hz) + 20 * math.log10(4 * math.pi / c)
    print(f"FSPL (dB) = 20*log10({d_m:.0f} m) + 20*log10({f_Hz:.0e} Hz) + 20*log10(4*pi/c)")
    print(f"FSPL = {fspl_dB:.2f} dB\n")

    print("--- Step 3: Calculate Received Signal Power (S) at Rx Antenna Terminals ---")
    S_dBm = eirp_dBm - fspl_dB + G_ant_rx_dB
    print(f"Signal Power (dBm) = EIRP - FSPL + Rx Antenna Gain")
    print(f"S = {eirp_dBm:.2f} dBm - {fspl_dB:.2f} dB + {G_ant_rx_dB} dBi = {S_dBm:.2f} dBm\n")

    print("--- Step 4: Calculate Receiver System Noise Figure (NF) ---")
    # Using Friis formula for cascaded noise figure.
    # The 'receiver' starts after the antenna gain is applied. First component is the antenna loss.
    # Stage 1: Rx Antenna Loss
    F1 = db_to_linear(L_ant_rx_dB)
    G1 = 1 / F1
    # Stage 2: Input Filter
    F2 = db_to_linear(L_filter_in_dB)
    G2 = 1 / F2
    # Stage 3: LNA
    F3 = db_to_linear(LNA_NF_dB)
    G3 = db_to_linear(LNA_gain_dB)
    # Stage 4: Mixer (NF = Loss, Gain = -Loss)
    F4 = db_to_linear(Mixer_loss_dB)
    G4 = 1 / F4

    # Calculate cumulative gain at each stage's input
    G_cum_1 = 1
    G_cum_2 = G1
    G_cum_3 = G1 * G2
    G_cum_4 = G1 * G2 * G3
    
    # Friis Formula: F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
    F_rx = F1 + (F2 - 1) / G_cum_2 + (F3 - 1) / G_cum_3 + (F4 - 1) / G_cum_4
    # Subsequent stages are ignored as their noise contribution is negligible due to the LNA's high gain.
    
    NF_rx_dB = linear_to_db(F_rx)
    print(f"The receiver noise figure is calculated by cascading the noise figures of its components.")
    print(f"Receiver System Noise Figure (NF) = {NF_rx_dB:.2f} dB\n")

    print("--- Step 5: Calculate Total Noise Power (N) at Receiver Input ---")
    # Noise power in dBm = Thermal Noise Density + 10*log10(Bandwidth) + Noise Figure
    noise_density_dBm_per_Hz = linear_to_db(k * T * 1000)
    N_dBm = noise_density_dBm_per_Hz + 10 * math.log10(B_Hz) + NF_rx_dB
    print(f"Noise Power (dBm) = Thermal Noise Density + 10*log10(Bandwidth) + NF")
    print(f"N = {noise_density_dBm_per_Hz:.2f} dBm/Hz + {10*math.log10(B_Hz):.2f} dB-Hz + {NF_rx_dB:.2f} dB = {N_dBm:.2f} dBm\n")

    print("--- Step 6: Calculate Final Signal-to-Noise Ratio (SNR) ---")
    snr_dB = S_dBm - N_dBm
    print("The final SNR is the difference between the received signal power and the total noise power.")
    print(f"SNR (dB) = Signal Power (dBm) - Noise Power (dBm)")
    print(f"SNR (dB) = {S_dBm:.2f} dBm - ({N_dBm:.2f} dBm)")
    print(f"Final SNR = {snr_dB:.2f} dB")
    
    return snr_dB

if __name__ == '__main__':
    final_snr = calculate_snr()
    # The final answer is wrapped according to the instruction.
    print(f"\n<<<{final_snr:.2f}>>>")