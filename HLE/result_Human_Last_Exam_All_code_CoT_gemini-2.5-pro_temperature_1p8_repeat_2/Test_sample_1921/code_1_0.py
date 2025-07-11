import math

def main():
    """
    Calculates the SNR of a microwave link based on provided system parameters.
    """
    # Given Parameters from the problem description
    tx_power_dbm = 30.0
    tx_ant_gain_db = 20.0
    tx_total_losses_db = 1.0 + 1.0 + 1.0  # Antenna Loss + Filter Loss + Cable Loss

    freq_ghz = 24.0
    distance_km = 10.0
    signal_bw_khz = 100.0

    rx_ant_gain_db = 1.0
    rx_pre_lna_losses_db = 0.5 + 1.0  # Antenna Loss + Input Filter Loss
    
    lna_gain_db = 36.0
    lna_nf_db = 2.0
    
    mixer_nf_db = 9.0  # Mixer conversion loss also defines its noise figure
    
    temp_k = 300.0
    k_boltzmann = 1.38e-23  # Boltzmann's constant in J/K

    # Step 1: Calculate Transmitter EIRP
    eirp_dbm = tx_power_dbm + tx_ant_gain_db - tx_total_losses_db

    # Step 2: Calculate Free Space Path Loss (FSPL)
    freq_hz = freq_ghz * 1e9
    distance_m = distance_km * 1e3
    # Using the FSPL formula: 20*log10(d_m) + 20*log10(f_hz) - 147.55
    fspl_db = 20 * math.log10(distance_m) + 20 * math.log10(freq_hz) - 147.55

    # Step 3: Calculate Received Signal Power (S_in) at the receiver input
    s_in_dbm = eirp_dbm - fspl_db + rx_ant_gain_db

    # Step 4: Calculate Receiver System Noise Figure (NF_sys)
    # Convert necessary values from dB to linear scale for Friis formula
    lna_gain_lin = 10**(lna_gain_db / 10)
    lna_nf_lin = 10**(lna_nf_db / 10)
    mixer_nf_lin = 10**(mixer_nf_db / 10)
    
    # Friis formula for cascaded noise figure of the active chain (LNA + Mixer)
    # F_active = F_lna + (F_mixer - 1) / G_lna
    # Contributions from later stages are negligible due to high LNA gain.
    f_active_chain_lin = lna_nf_lin + (mixer_nf_lin - 1) / lna_gain_lin
    nf_active_chain_db = 10 * math.log10(f_active_chain_lin)
    
    # Total system noise figure is the sum of pre-LNA losses and the active chain's NF
    nf_sys_db = rx_pre_lna_losses_db + nf_active_chain_db

    # Step 5: Calculate Input-Referred Noise Power (N_in)
    signal_bw_hz = signal_bw_khz * 1e3
    # Thermal noise power in dBm = 10*log10(k * T * B * 1000 mW/W)
    thermal_noise_power_dbm = 10 * math.log10(k_boltzmann * temp_k * signal_bw_hz * 1000)
    # Input-Referred Noise is the thermal noise plus the noise added by the receiver system
    n_in_dbm = thermal_noise_power_dbm + nf_sys_db

    # Step 6: Calculate Final SNR
    snr_db = s_in_dbm - n_in_dbm

    # Print the final equation with all numbers
    print("--- Final SNR Calculation ---")
    print(f"The SNR is calculated as the Input Signal Power minus the Input-Referred Noise Power.")
    print(f"SNR (dB) = S_in (dBm) - N_in (dBm)")
    print(f"Final Equation: {snr_db:.2f} dB = {s_in_dbm:.2f} dBm - ({n_in_dbm:.2f} dBm)")
    print(f"\nThe resulting SNR of the received signal is {snr_db:.2f} dB.")

if __name__ == '__main__':
    main()
<<<28.27>>>