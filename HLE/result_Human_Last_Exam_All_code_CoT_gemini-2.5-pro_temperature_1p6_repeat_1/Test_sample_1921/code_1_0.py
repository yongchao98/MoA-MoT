import math

def db_to_linear(db):
    """Converts a value from dB to linear scale."""
    return 10**(db / 10)

def linear_to_db(linear):
    """Converts a value from linear to dB scale."""
    return 10 * math.log10(linear)

def main():
    # --- Given Parameters ---
    # Transmitter (Tx)
    tx_power_dbm = 30.0
    tx_ant_gain_db = 20.0
    tx_ant_loss_db = 1.0
    tx_filter_loss_db = 1.0
    tx_cable_loss_db = 1.0

    # Path
    freq_ghz = 24.0
    distance_km = 10.0

    # Receiver (Rx)
    rx_ant_gain_db = 1.0
    rx_ant_loss_db = 0.5
    rx_input_filter_loss_db = 1.0
    lna_gain_db = 36.0
    lna_nf_db = 2.0
    mixer_loss_db = 9.0
    if_filter_loss_db = 1.0
    if_amp_gain_db = 23.0
    if_amp_nf_db = 0.0  # Negligible NF
    output_filter_loss_db = 1.0

    # System
    bandwidth_hz = 100 * 1000  # 100 kHz
    temp_k = 300.0  # Ambient temperature in Kelvin
    BOLTZMANN_K = 1.380649e-23  # Boltzmann's constant

    print("--- Step 1: Calculate Transmitter EIRP ---")
    tx_losses_db = tx_ant_loss_db + tx_filter_loss_db + tx_cable_loss_db
    eirp_dbm = tx_power_dbm + tx_ant_gain_db - tx_losses_db
    print(f"Tx Power: {tx_power_dbm} dBm")
    print(f"Tx Antenna Gain: {tx_ant_gain_db} dB")
    print(f"Total Tx Losses: {tx_losses_db} dB")
    print(f"EIRP = {tx_power_dbm} + {tx_ant_gain_db} - {tx_losses_db} = {eirp_dbm:.2f} dBm\n")

    print("--- Step 2: Calculate Free Space Path Loss (FSPL) ---")
    freq_mhz = freq_ghz * 1000
    fspl_db = 20 * math.log10(distance_km) + 20 * math.log10(freq_mhz) + 32.45
    print(f"Frequency: {freq_ghz} GHz, Distance: {distance_km} km")
    print(f"FSPL = 20*log10({distance_km}) + 20*log10({freq_mhz}) + 32.45 = {fspl_db:.2f} dB\n")

    print("--- Step 3: Calculate Received Signal Power (S) ---")
    # This is the power at the input of the first electronic component (input filter)
    received_signal_dbm = eirp_dbm - fspl_db + rx_ant_gain_db - rx_ant_loss_db
    print("Signal power at receiver input is calculated as: EIRP - FSPL + Rx Ant Gain - Rx Ant Loss")
    print(f"S = {eirp_dbm:.2f} dBm - {fspl_db:.2f} dB + {rx_ant_gain_db} dB - {rx_ant_loss_db} dB")
    print(f"Signal Power (S_in): {received_signal_dbm:.2f} dBm\n")
    
    print("--- Step 4: Calculate Receiver Noise Temperature (T_e) ---")
    # Temperature of components is referenced to the ambient temperature given (300K)
    # T_e = T_comp * (F-1) for active components, T_e = T_phys * (L-1) for passive
    
    # Component values in linear scale and their noise temperatures
    # Stage 1: Input Filter
    l1_lin = db_to_linear(rx_input_filter_loss_db)
    t1_k = temp_k * (l1_lin - 1)
    g1_lin = 1 / l1_lin
    # Stage 2: LNA
    f2_lin = db_to_linear(lna_nf_db)
    t2_k = temp_k * (f2_lin - 1) # Using ambient temp as reference per problem statement
    g2_lin = db_to_linear(lna_gain_db)
    # Stage 3: Mixer
    l3_lin = db_to_linear(mixer_loss_db)
    t3_k = temp_k * (l3_lin - 1)
    g3_lin = 1 / l3_lin
    # Stage 4: IF Filter
    l4_lin = db_to_linear(if_filter_loss_db)
    t4_k = temp_k * (l4_lin - 1)
    g4_lin = 1 / l4_lin
    # Stage 5: IF Amp
    f5_lin = db_to_linear(if_amp_nf_db)
    t5_k = temp_k * (f5_lin - 1)
    g5_lin = db_to_linear(if_amp_gain_db)
    # Stage 6: Output Filter
    l6_lin = db_to_linear(output_filter_loss_db)
    t6_k = temp_k * (l6_lin - 1)

    # Cascaded Noise Temperature Calculation
    t_e_receiver = t1_k \
                + t2_k / g1_lin \
                + t3_k / (g1_lin * g2_lin) \
                + t4_k / (g1_lin * g2_lin * g3_lin) \
                + t5_k / (g1_lin * g2_lin * g3_lin * g4_lin) \
                + t6_k / (g1_lin * g2_lin * g3_lin * g4_lin * g5_lin)
    print(f"Receiver effective noise temperature (T_e) = {t_e_receiver:.2f} K\n")

    print("--- Step 5: Calculate Total Input Noise Power (N) ---")
    t_ant_k = temp_k # Antenna noise temperature assumed to be ambient
    t_sys_k = t_ant_k + t_e_receiver
    print(f"System Noise Temperature (T_sys) = T_antenna + T_e = {t_ant_k:.2f} K + {t_e_receiver:.2f} K = {t_sys_k:.2f} K")
    
    noise_power_watts = BOLTZMANN_K * t_sys_k * bandwidth_hz
    noise_power_dbm = linear_to_db(noise_power_watts * 1000)
    print("Noise power is calculated as: k * T_sys * Bandwidth")
    print(f"N = {BOLTZMANN_K:.4e} * {t_sys_k:.2f} * {bandwidth_hz}")
    print(f"Noise Power (N_in): {noise_power_dbm:.2f} dBm\n")

    print("--- Step 6: Calculate Final SNR ---")
    snr_db = received_signal_dbm - noise_power_dbm
    print("SNR is the difference between the signal power and the noise power in dB.")
    print(f"SNR (dB) = S (dBm) - N (dBm)")
    print(f"SNR (dB) = {received_signal_dbm:.2f} dBm - ({noise_power_dbm:.2f} dBm)")
    print(f"Final SNR = {snr_db:.2f} dB")
    
    return snr_db

if __name__ == "__main__":
    final_snr = main()
    print(f"\n<<<{final_snr:.2f}>>>")