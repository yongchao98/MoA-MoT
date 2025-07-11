import math

# A helper function for dB conversions
def from_db(x):
    """Converts from dB to linear scale."""
    return 10**(x / 10)

def main():
    # --- System Parameters ---
    P_tx_dBm = 30.0         # Transmitter Power in dBm
    G_tx_dB = 20.0          # Tx Antenna Gain in dB
    L_cable_tx_dB = 1.0     # Tx Cable loss in dB
    L_filt_tx_dB = 1.0      # Tx Output filter loss in dB
    L_ant_tx_dB = 1.0       # Tx Antenna loss in dB
    
    freq_Hz = 24e9          # Frequency in Hz (24 GHz)
    dist_m = 10e3           # Distance in meters (10 km)
    
    G_rx_dB = 1.0           # Rx Antenna gain in dB
    L_ant_rx_dB = 0.5       # Rx Antenna loss in dB
    
    # Receiver Chain components (in order)
    # 1. Input Filter
    L_filt_rx_in_dB = 1.0
    # 2. LNA
    G_lna_dB = 36.0
    NF_lna_dB = 2.0
    # 3. Mixer
    L_mix_dB = 9.0
    NF_mix_dB = 9.0 # For passive mixers, NF is assumed equal to conversion loss
    
    # System Constants
    BW_signal_Hz = 100e3    # Signal Bandwidth in Hz (100 kHz)
    T_kelvin = 300.0        # Ambient Temperature in Kelvin
    k_boltzmann = 1.38e-23  # Boltzmann's constant in J/K
    c_light = 3e8           # Speed of light in m/s

    # --- Step 1: Calculate Effective Isotropic Radiated Power (EIRP) ---
    print("--- Calculating Signal Power ---")
    EIRP_dBm = P_tx_dBm - L_cable_tx_dB - L_filt_tx_dB + G_tx_dB - L_ant_tx_dB
    print(f"1. Effective Isotropic Radiated Power (EIRP):")
    print(f"   EIRP = P_tx - L_cable - L_filter + G_tx - L_ant_tx")
    print(f"   EIRP = {P_tx_dBm} dBm - {L_cable_tx_dB} dB - {L_filt_tx_dB} dB + {G_tx_dB} dB - {L_ant_tx_dB} dB = {EIRP_dBm:.2f} dBm\n")

    # --- Step 2: Calculate Free Space Path Loss (FSPL) ---
    # Using the formula FSPL (dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
    fspl_dB = 20 * math.log10(dist_m) + 20 * math.log10(freq_Hz) + 20 * math.log10(4 * math.pi / c_light)
    print(f"2. Free Space Path Loss (FSPL) at {freq_Hz/1e9} GHz over {dist_m/1e3} km:")
    print(f"   FSPL = {fspl_dB:.2f} dB\n")

    # --- Step 3: Calculate Received Signal Power (S) ---
    # This is the power at the input terminals of the receiver electronics (after Rx antenna but before filters)
    S_rx_dBm = EIRP_dBm - fspl_dB + G_rx_dB - L_ant_rx_dB
    print(f"3. Received Signal Power (S) at receiver input:")
    print(f"   S = EIRP - FSPL + G_rx - L_ant_rx")
    print(f"   S = {EIRP_dBm:.2f} dBm - {fspl_dB:.2f} dB + {G_rx_dB} dB - {L_ant_rx_dB} dB = {S_rx_dBm:.2f} dBm\n")

    # --- Step 4: Calculate System Noise Figure (NF_sys) ---
    print("--- Calculating Noise Power ---")
    # Reference point for NF is the same as for signal S: receiver input.
    # The chain for NF calculation starts with the first component, the input filter.
    # We use the Friis formula: F_sys = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
    F1_lin = from_db(L_filt_rx_in_dB)
    G1_lin = 1 / F1_lin

    F2_lin = from_db(NF_lna_dB)
    G2_lin = from_db(G_lna_dB)
    
    F3_lin = from_db(NF_mix_dB)
    # Contribution from later stages is negligible due to high LNA gain.

    F_sys_lin = F1_lin + (F2_lin - 1) / G1_lin + (F3_lin - 1) / (G1_lin * G2_lin)
    NF_sys_dB = 10 * math.log10(F_sys_lin)
    print(f"4. Receiver System Noise Figure (NF_sys):")
    print(f"   Calculated using Friis formula for the receiver chain (filter, LNA, mixer).")
    print(f"   NF_sys = {NF_sys_dB:.2f} dB\n")

    # --- Step 5: Calculate Total Noise Power (N) ---
    # First, calculate thermal noise in the signal bandwidth
    N_thermal_watts = k_boltzmann * T_kelvin * BW_signal_Hz
    N_thermal_dBm = 10 * math.log10(N_thermal_watts / 0.001)
    
    # Total noise is thermal noise plus noise from the receiver electronics (NF)
    N_total_dBm = N_thermal_dBm + NF_sys_dB
    print(f"5. Total Noise Power (N) in {BW_signal_Hz/1e3} kHz bandwidth:")
    print(f"   Thermal Noise = 10*log10(k*T*B / 1mW) = {N_thermal_dBm:.2f} dBm")
    print(f"   Total Noise (N) = Thermal Noise + NF_sys")
    print(f"   N = {N_thermal_dBm:.2f} dBm + {NF_sys_dB:.2f} dB = {N_total_dBm:.2f} dBm\n")

    # --- Step 6: Calculate Final SNR ---
    SNR_dB = S_rx_dBm - N_total_dBm
    print("--- Calculating Final SNR ---")
    print(f"6. Final Signal-to-Noise Ratio (SNR):")
    print(f"   SNR = Signal Power (S) - Noise Power (N)")
    print(f"   SNR = {S_rx_dBm:.2f} dBm - ({N_total_dBm:.2f} dBm)")
    print(f"   Resulting SNR = {SNR_dB:.2f} dB")
    
    return SNR_dB

if __name__ == '__main__':
    final_snr = main()
    # The final answer is enclosed in <<<>>>
    print(f"\n<<<Result: {final_snr:.2f} dB>>>")

<<<28.27>>>