import math

def solve_snr():
    """
    Calculates the Signal-to-Noise Ratio (SNR) for a microwave link based on given parameters.
    """

    # --- Step 0: Define constants and initial parameters ---
    # Given parameters
    Tx_power_dBm = 30.0
    Tx_ant_gain_dB = 20.0
    Tx_ant_loss_dB = 1.0
    Tx_filter_loss_dB = 1.0
    Tx_cable_loss_dB = 1.0
    freq_Hz = 24e9  # 24 GHz
    distance_m = 10e3  # 10 km
    bandwidth_Hz = 100e3  # 100 kHz
    Rx_ant_gain_dB = 1.0
    Rx_ant_loss_dB = 0.5
    Rx_filter_loss_dB = 1.0
    LNA_gain_dB = 36.0
    LNA_nf_dB = 2.0
    Mixer_loss_dB = 9.0
    Mixer_nf_dB = 9.0 # Assuming mixer noise figure equals its conversion loss
    IF_filter_loss_dB = 1.0
    IF_filter_nf_dB = 1.0 # Assuming passive filter NF equals its loss
    IF_amp_gain_dB = 23.0
    IF_amp_nf_dB = 0.0 # Negligible NF
    
    # Physical constants
    k_boltzmann = 1.38e-23  # J/K
    temp_K = 300.0  # Kelvin
    c_light = 3.0e8 # m/s

    print("--- Link Budget Calculation ---")
    
    # --- Step 1: Calculate Transmitter EIRP ---
    tx_total_loss_dB = Tx_ant_loss_dB + Tx_filter_loss_dB + Tx_cable_loss_dB
    eirp_dBm = Tx_power_dBm + Tx_ant_gain_dB - tx_total_loss_dB
    print(f"1. Transmitter EIRP:")
    print(f"   EIRP (dBm) = Tx Power ({Tx_power_dBm} dBm) + Ant Gain ({Tx_ant_gain_dB} dB) - Tx Losses ({tx_total_loss_dB} dB) = {eirp_dBm:.2f} dBm")
    
    # --- Step 2: Calculate Free Space Path Loss (FSPL) ---
    # FSPL (dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
    fspl_dB = 20 * math.log10(distance_m) + 20 * math.log10(freq_Hz) + 20 * math.log10(4 * math.pi / c_light)
    print(f"\n2. Free Space Path Loss (FSPL):")
    print(f"   FSPL for {distance_m/1000} km at {freq_Hz/1e9} GHz = {fspl_dB:.2f} dB")
    
    # --- Step 3: Calculate Received Signal Power (S) ---
    signal_power_dBm = eirp_dBm - fspl_dB + Rx_ant_gain_dB
    print(f"\n3. Received Signal Power (S) at antenna port:")
    print(f"   S (dBm) = EIRP ({eirp_dBm:.2f} dBm) - FSPL ({fspl_dB:.2f} dB) + Rx Ant Gain ({Rx_ant_gain_dB} dB) = {signal_power_dBm:.2f} dBm")

    # --- Step 4: Calculate Total System Noise Figure (NF) ---
    print(f"\n4. Receiver Cascaded Noise Figure (NF):")
    
    # Helper to convert dB to linear scale
    def db_to_linear(db_val):
        return 10**(db_val / 10.0)

    # Receiver chain components [name, gain_dB, noise_figure_dB]
    # For passive components, Gain = -Loss and NF = Loss
    rx_chain = [
        ("Rx Ant Loss", -Rx_ant_loss_dB, Rx_ant_loss_dB),
        ("Rx Filter", -Rx_filter_loss_dB, Rx_filter_loss_dB),
        ("LNA", LNA_gain_dB, LNA_nf_dB),
        ("Mixer", -Mixer_loss_dB, Mixer_nf_dB),
        ("IF Filter", -IF_filter_loss_dB, IF_filter_nf_dB),
        ("IF Amp", IF_amp_gain_dB, IF_amp_nf_dB)
    ]

    total_nf_lin = 0.0
    cascaded_gain_lin = 1.0

    # Apply Friis formula for cascaded noise figure
    for i, (name, gain_db, nf_db) in enumerate(rx_chain):
        gain_lin = db_to_linear(gain_db)
        nf_lin = db_to_linear(nf_db)
        
        if i == 0:
            stage_contribution = nf_lin
        else:
            stage_contribution = (nf_lin - 1) / cascaded_gain_lin
        
        total_nf_lin += stage_contribution
        cascaded_gain_lin *= gain_lin
        # print(f"   - After {name}: Stage Contribution={stage_contribution:.4f}, Cumulative NF={10*math.log10(total_nf_lin):.2f} dB")
    
    total_nf_dB = 10 * math.log10(total_nf_lin)
    print(f"   Total System Noise Figure = {total_nf_dB:.2f} dB")

    # --- Step 5: Calculate Total Input-Referred Noise Power (N) ---
    thermal_noise_watts = k_boltzmann * temp_K * bandwidth_Hz
    thermal_noise_dBm = 10 * math.log10(thermal_noise_watts / 0.001)
    
    total_noise_power_dBm = thermal_noise_dBm + total_nf_dB
    print(f"\n5. Total Input-Referred Noise Power (N):")
    print(f"   Thermal Noise in {bandwidth_Hz/1000} kHz BW ({thermal_noise_dBm:.2f} dBm) + System NF ({total_nf_dB:.2f} dB) = {total_noise_power_dBm:.2f} dBm")

    # --- Step 6: Calculate Final SNR ---
    final_snr_dB = signal_power_dBm - total_noise_power_dBm
    print(f"\n--- Final Result ---")
    print(f"6. Resulting Signal-to-Noise Ratio (SNR):")
    print(f"   SNR (dB) = Signal Power (S) - Total Noise Power (N)")
    print(f"   SNR (dB) = {signal_power_dBm:.2f} dBm - ({total_noise_power_dBm:.2f} dBm)")
    print(f"   SNR (dB) = {final_snr_dB:.2f} dB")
    
    return final_snr_dB

# Execute the function and capture the final answer
final_answer = solve_snr()

# Print the final answer in the required format
# print(f"\n<<<{final_answer:.2f}>>>")
print(f"<<<{round(final_answer, 2)}>>>")
