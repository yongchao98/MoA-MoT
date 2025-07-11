import math

def db_to_linear(db_val):
    """Converts a value from dB to linear scale."""
    return 10**(db_val / 10)

def linear_to_db(linear_val):
    """Converts a value from linear to dB scale."""
    return 10 * math.log10(linear_val)

def main():
    # --- Input Parameters ---
    # Transmitter
    tx_power_dbm = 30.0
    tx_ant_gain_db = 20.0
    tx_losses_db = 1.0 + 1.0 + 1.0  # ant_loss + filter_loss + cable_loss

    # Path
    frequency_ghz = 24.0
    distance_km = 10.0

    # Receiver
    rx_ant_gain_db = 1.0
    # Receiver components in order [gain_dB, noise_figure_dB]
    # For passive components, Gain = -Loss and NF = Loss
    rx_chain = [
        [-0.5, 0.5],  # Rx Antenna Loss
        [-1.0, 1.0],  # Rx Input Filter
        [36.0, 2.0],  # LNA
        [-9.0, 9.0],  # Mixer (Conversion Loss and NF)
        [-1.0, 1.0],  # IF Filter
        [23.0, 0.0],  # IF Amplifier (negligible NF)
        [-1.0, 1.0]   # Output Filter
    ]

    # System
    temp_k = 300.0
    bandwidth_hz = 100e3
    k_boltzmann = 1.380649e-23  # J/K

    print("--- Link Budget Calculation ---")

    # 1. Calculate EIRP (Effective Isotropic Radiated Power)
    eirp_dbm = tx_power_dbm + tx_ant_gain_db - tx_losses_db
    print(f"1. Tx Power: {tx_power_dbm} dBm, Tx Antenna Gain: {tx_ant_gain_db} dB, Tx Losses: {tx_losses_db} dB")
    print(f"   => EIRP = {tx_power_dbm} + {tx_ant_gain_db} - {tx_losses_db} = {eirp_dbm:.2f} dBm\n")

    # 2. Calculate FSPL (Free Space Path Loss)
    fspl_db = 92.45 + 20 * math.log10(frequency_ghz) + 20 * math.log10(distance_km)
    print(f"2. Frequency: {frequency_ghz} GHz, Distance: {distance_km} km")
    print(f"   => Free Space Path Loss (FSPL) = {fspl_db:.2f} dB\n")

    # 3. Calculate Received Signal Power (at Rx antenna input)
    received_signal_dbm = eirp_dbm - fspl_db + rx_ant_gain_db
    print(f"3. EIRP: {eirp_dbm:.2f} dBm, FSPL: {fspl_db:.2f} dB, Rx Antenna Gain: {rx_ant_gain_db} dB")
    print(f"   => Received Signal Power = {eirp_dbm:.2f} - {fspl_db:.2f} + {rx_ant_gain_db} = {received_signal_dbm:.2f} dBm\n")

    # 4. Calculate Receiver System Noise Figure (using Friis formula)
    total_noise_factor_lin = 0
    cumulative_gain_lin = 1.0
    
    components_lin = []
    for gain_db, nf_db in rx_chain:
        components_lin.append({'gain': db_to_linear(gain_db), 'nf': db_to_linear(nf_db)})

    # Friis formula: F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
    f1 = components_lin[0]['nf']
    g1 = components_lin[0]['gain']
    f_sys_lin = f1
    gain_cascade_lin = g1
    
    for i in range(1, len(components_lin)):
        fi = components_lin[i]['nf']
        gi = components_lin[i]['gain']
        f_sys_lin += (fi - 1) / gain_cascade_lin
        gain_cascade_lin *= gi

    system_nf_db = linear_to_db(f_sys_lin)
    print(f"4. Calculating cascaded noise figure for the receiver chain...")
    print(f"   => Total Receiver System Noise Figure = {system_nf_db:.2f} dB\n")

    # 5. Calculate Total Input-Referred Noise Power
    thermal_noise_watts = k_boltzmann * temp_k * bandwidth_hz
    thermal_noise_dbm = linear_to_db(thermal_noise_watts * 1000)
    print(f"5. Temperature: {temp_k} K, Bandwidth: {bandwidth_hz/1e3} kHz")
    print(f"   => Thermal Noise Floor = {thermal_noise_dbm:.2f} dBm")
    
    total_noise_dbm = thermal_noise_dbm + system_nf_db
    print(f"   Total Input Noise Power = Thermal Noise + System NF")
    print(f"   => Total Input Noise Power = {thermal_noise_dbm:.2f} dBm + {system_nf_db:.2f} dB = {total_noise_dbm:.2f} dBm\n")
    
    # 6. Calculate Final SNR
    snr_db = received_signal_dbm - total_noise_dbm
    print("--- Final Result ---")
    print(f"6. SNR (dB) = Received Signal Power (dBm) - Total Input Noise Power (dBm)")
    print(f"   Final SNR = {received_signal_dbm:.2f} - ({total_noise_dbm:.2f})")
    print(f"   Resulting SNR = {snr_db:.2f} dB")


if __name__ == "__main__":
    main()