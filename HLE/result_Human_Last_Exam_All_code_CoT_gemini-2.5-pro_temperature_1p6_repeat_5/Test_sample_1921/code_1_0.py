import math

# --- Helper functions ---
def to_linear(db_value):
    """Converts a dB value to a linear ratio."""
    return 10**(db_value / 10)

def to_db(linear_value):
    """Converts a linear ratio to dB."""
    return 10 * math.log10(linear_value)

def main():
    """
    Calculates the Signal-to-Noise Ratio (SNR) for the given microwave link.
    """
    # --- Step 0: Define constants and given parameters ---

    # Constants
    BOLTZMANN_K = 1.38e-23  # Boltzmann's constant in J/K
    TEMP_K = 300           # Ambient temperature in Kelvin
    
    # Link Parameters
    FREQ_GHZ = 24.0        # Operating frequency in GHz
    DISTANCE_KM = 10.0     # Distance between Tx and Rx in km
    BANDWIDTH_HZ = 100e3   # Signal bandwidth in Hz

    # Transmitter (Tx) Parameters
    TX_POWER_DBM = 30.0
    TX_ANTENNA_GAIN_DB = 20.0
    TX_ANTENNA_LOSS_DB = 1.0
    TX_FILTER_LOSS_DB = 1.0
    TX_CABLE_LOSS_DB = 1.0

    # Receiver (Rx) Parameters
    RX_ANTENNA_GAIN_DB = 1.0
    RX_ANTENNA_LOSS_DB = 0.5
    
    # Receiver chain components [(gain_db, noise_figure_db)]
    # For passive components (filters, lossy mixers), Noise Figure = Loss
    rx_chain = [
        (-1.0, 1.0),   # Rx Input Filter (1dB loss, 1dB NF)
        (36.0, 2.0),   # LNA (36dB gain, 2dB NF)
        (-9.0, 9.0),   # Mixer (9dB conversion loss, 9dB NF)
        (-1.0, 1.0),   # IF Filter (1dB loss, 1dB NF)
        (23.0, 0.0),   # IF Amplifier (23dB gain, negligible NF)
        (-1.0, 1.0)    # Output Filter (1dB loss, 1dB NF)
    ]

    # --- Step 1: Calculate Transmitter EIRP (Effective Isotropic Radiated Power) ---
    tx_total_loss_db = TX_ANTENNA_LOSS_DB + TX_FILTER_LOSS_DB + TX_CABLE_LOSS_DB
    eirp_dbm = TX_POWER_DBM + TX_ANTENNA_GAIN_DB - tx_total_loss_db

    # --- Step 2: Calculate Free Space Path Loss (FSPL) ---
    # Formula for FSPL in dB with frequency in GHz and distance in km:
    # FSPL(dB) = 20*log10(distance_km) + 20*log10(frequency_GHz) + 92.45
    fspl_db = 20 * math.log10(DISTANCE_KM) + 20 * math.log10(FREQ_GHZ) + 92.45

    # --- Step 3: Calculate Received Signal Power at Rx input ---
    # This is the power at the input to the first electronic component, after the Rx antenna.
    p_rx_input_dbm = eirp_dbm - fspl_db + RX_ANTENNA_GAIN_DB - RX_ANTENNA_LOSS_DB

    # --- Step 4: Calculate Thermal Noise Power in the signal bandwidth ---
    noise_power_watts = BOLTZMANN_K * TEMP_K * BANDWIDTH_HZ
    # Convert noise power from Watts to dBm
    noise_power_dbm = 10 * math.log10(noise_power_watts / 0.001)

    # --- Step 5: Calculate the cascaded System Noise Figure (NF) ---
    # Use the Friis formula: F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
    total_noise_factor_lin = 0
    cascaded_gain_lin = 1.0
    
    # Unpack first stage to initialize the total noise factor
    g1_db, nf1_db = rx_chain[0]
    f1_lin = to_linear(nf1_db)
    g1_lin = to_linear(g1_db)
    total_noise_factor_lin = f1_lin
    cascaded_gain_lin = g1_lin

    # Loop through the rest of the stages
    for g_db, nf_db in rx_chain[1:]:
        f_lin = to_linear(nf_db)
        g_lin = to_linear(g_db)
        total_noise_factor_lin += (f_lin - 1) / cascaded_gain_lin
        cascaded_gain_lin *= g_lin
        
    system_nf_db = to_db(total_noise_factor_lin)

    # --- Step 6: Calculate the final SNR ---
    # The output SNR is the input SNR minus the system's Noise Figure.
    # SNR_out(dB) = SNR_in(dB) - NF_sys(dB)
    # SNR_in(dB) = P_rx_input(dBm) - N_thermal(dBm)
    final_snr_db = p_rx_input_dbm - noise_power_dbm - system_nf_db

    # --- Print the final result ---
    print("This script calculates the final Signal-to-Noise Ratio (SNR).")
    print("The calculation is: SNR(dB) = P_rx_input(dBm) - N_thermal(dBm) - NF_system(dB)\n")
    print(f"Calculated values:")
    print(f"  - Received Signal Power (P_rx_input): {p_rx_input_dbm:.2f} dBm")
    print(f"  - Thermal Noise Power (N_thermal): {noise_power_dbm:.2f} dBm")
    print(f"  - System Noise Figure (NF_system): {system_nf_db:.2f} dB\n")
    print("Final SNR Equation:")
    print(f"SNR = {p_rx_input_dbm:.2f} dBm - ({noise_power_dbm:.2f} dBm) - {system_nf_db:.2f} dB")
    print(f"Resulting SNR = {final_snr_db:.2f} dB")

if __name__ == "__main__":
    main()
<<<28.27>>>