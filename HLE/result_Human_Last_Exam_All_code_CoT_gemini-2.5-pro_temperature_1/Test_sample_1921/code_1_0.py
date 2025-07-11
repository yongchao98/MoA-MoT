import math

def calculate_snr():
    """
    Calculates the Signal-to-Noise Ratio (SNR) for the given microwave link budget.
    """
    # --- Given values ---
    Tx_power_dBm = 30.0
    Tx_gain_dB = 20.0
    Tx_ant_loss_dB = 1.0
    Tx_filter_loss_dB = 1.0
    Tx_cable_loss_dB = 1.0

    freq_GHz = 24.0
    dist_km = 10.0

    Rx_gain_dB = 1.0
    Rx_ant_loss_dB = 0.5
    Rx_filter1_loss_dB = 1.0
    LNA_gain_dB = 36.0
    LNA_NF_dB = 2.0
    
    BW_Hz = 100 * 1000  # 100 kHz Signal Bandwidth

    print("This script calculates the final SNR of the received signal.\n")
    
    # --- 1. Calculate EIRP (Effective Isotropic Radiated Power) ---
    EIRP_dBm = Tx_power_dBm - Tx_filter_loss_dB - Tx_cable_loss_dB - Tx_ant_loss_dB + Tx_gain_dB
    print("Step 1: Calculate Effective Isotropic Radiated Power (EIRP)")
    print(f"EIRP (dBm) = Tx Power ({Tx_power_dBm} dBm) - Tx Filter Loss ({Tx_filter_loss_dB} dB) - Tx Cable Loss ({Tx_cable_loss_dB} dB) - Tx Antenna Loss ({Tx_ant_loss_dB} dB) + Tx Antenna Gain ({Tx_gain_dB} dB)")
    print(f"Result: {Tx_power_dBm} - {Tx_filter_loss_dB} - {Tx_cable_loss_dB} - {Tx_ant_loss_dB} + {Tx_gain_dB} = {EIRP_dBm:.2f} dBm\n")

    # --- 2. Calculate Free Space Path Loss (FSPL) ---
    # Using the formula: FSPL (dB) = 32.45 + 20*log10(f_MHz) + 20*log10(d_km)
    freq_MHz = freq_GHz * 1000
    FSPL_dB = 32.45 + 20 * math.log10(freq_MHz) + 20 * math.log10(dist_km)
    print("Step 2: Calculate Free Space Path Loss (FSPL)")
    print(f"FSPL (dB) = 32.45 + 20*log10(Frequency in MHz) + 20*log10(Distance in km)")
    print(f"Result: 32.45 + 20*log10({freq_MHz}) + 20*log10({dist_km}) = {FSPL_dB:.2f} dB\n")

    # --- 3. Calculate Received Signal Power (S_in) at Rx Antenna Input ---
    S_in_dBm = EIRP_dBm - FSPL_dB + Rx_gain_dB
    print("Step 3: Calculate Received Signal Power (S_in) at Rx Antenna Input")
    print(f"S_in (dBm) = EIRP ({EIRP_dBm:.2f} dBm) - FSPL ({FSPL_dB:.2f} dB) + Rx Antenna Gain ({Rx_gain_dB} dB)")
    print(f"Result: {EIRP_dBm:.2f} - {FSPL_dB:.2f} + {Rx_gain_dB} = {S_in_dBm:.2f} dBm\n")

    # --- 4. Calculate System Noise Figure (NF_sys) ---
    # NF_sys is approximated by the losses before the LNA plus the LNA's own NF.
    loss_before_LNA_dB = Rx_ant_loss_dB + Rx_filter1_loss_dB
    NF_sys_dB = loss_before_LNA_dB + LNA_NF_dB
    print("Step 4: Calculate System Noise Figure (NF_sys)")
    print(f"NF_sys (dB) ≈ Loss before LNA + LNA NF")
    print(f"Result: ({Rx_ant_loss_dB} + {Rx_filter1_loss_dB}) + {LNA_NF_dB} = {NF_sys_dB:.2f} dB\n")

    # --- 5. Calculate Input-Referred Noise Power (N_in) ---
    # Thermal noise power in dBm = -174 dBm/Hz + 10*log10(Bandwidth in Hz)
    thermal_noise_dBm = -174 + 10 * math.log10(BW_Hz)
    N_in_dBm = thermal_noise_dBm + NF_sys_dB
    print("Step 5: Calculate Input-Referred Noise Power (N_in)")
    print(f"N_in (dBm) = (Thermal Noise @ 300K) + System NF")
    print(f"Thermal Noise Part = -174 dBm/Hz + 10*log10({BW_Hz} Hz) = {thermal_noise_dBm:.2f} dBm")
    print(f"Result: {thermal_noise_dBm:.2f} dBm + {NF_sys_dB:.2f} dB = {N_in_dBm:.2f} dBm\n")

    # --- 6. Calculate Final SNR ---
    SNR_dB = S_in_dBm - N_in_dBm
    print("Step 6: Calculate Final Signal-to-Noise Ratio (SNR)")
    print(f"SNR (dB) = Input Signal Power (S_in) - Input Noise Power (N_in)")
    print(f"Final Equation: {S_in_dBm:.2f} dBm - ({N_in_dBm:.2f} dBm) = {SNR_dB:.2f} dB")
    
    print(f"\nFinal Answer: The resulting SNR of the received signal is {SNR_dB:.2f} dB.")

if __name__ == '__main__':
    calculate_snr()
    print("<<<28.45>>>")
