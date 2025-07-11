import math

def db_to_linear(db):
    """Converts a dB value to a linear scale."""
    return 10**(db / 10)

def linear_to_db(lin):
    """Converts a linear value to a dB scale."""
    return 10 * math.log10(lin)

def main():
    # --- System Parameters ---
    # Constants
    k = 1.380649e-23  # Boltzmann's constant in J/K
    T = 300            # Ambient temperature in Kelvin
    c = 299792458      # Speed of light in m/s

    # Transmitter (Tx)
    Ptx_dBm = 30.0       # Transmitter power
    Gtx_dB = 20.0        # Tx antenna gain
    Ltx_ant_dB = 1.0     # Tx antenna loss
    Ltx_filt_dB = 1.0    # Tx filter loss
    Ltx_cable_dB = 1.0   # Tx cable loss
    freq_Hz = 24e9       # Operating frequency (24 GHz)
    B_Hz = 100e3         # Modulated signal bandwidth (100 kHz)

    # Path
    dist_m = 10e3        # Distance (10 km)

    # Receiver (Rx)
    Grx_dB = 1.0         # Rx antenna gain
    
    # Receiver chain components [Gain (dB), Noise Figure (dB)]
    # For passive components, Loss = Noise Figure, and Gain = -Loss.
    rx_chain = [
        {'name': 'Rx Antenna',      'gain_db': -0.5, 'nf_db': 0.5},
        {'name': 'Rx Input Filter', 'gain_db': -1.0, 'nf_db': 1.0},
        {'name': 'LNA',             'gain_db': 36.0, 'nf_db': 2.0},
        {'name': 'Mixer',           'gain_db': -9.0, 'nf_db': 9.0}, # Conversion Loss acts as Gain and NF
        {'name': 'IF Filter',       'gain_db': -1.0, 'nf_db': 1.0},
        {'name': 'IF Amplifier',    'gain_db': 23.0, 'nf_db': 0.0}, # Negligible NF
        {'name': 'Output Filter',   'gain_db': -1.0, 'nf_db': 1.0}
    ]

    print("--- Step 1: Calculate Received Signal Power (Prx) ---")
    
    # Calculate total transmitter losses
    Ltx_total_dB = Ltx_ant_dB + Ltx_filt_dB + Ltx_cable_dB
    print(f"Total Transmitter Loss (Ltx) = {Ltx_ant_dB} dB + {Ltx_filt_dB} dB + {Ltx_cable_dB} dB = {Ltx_total_dB:.2f} dB")
    
    # Calculate Free Space Path Loss (Lfs)
    # Lfs (dB) = 20 * log10(d) + 20 * log10(f) + 20 * log10(4*pi/c)
    lfs_val = (4 * math.pi * dist_m * freq_Hz) / c
    Lfs_dB = 20 * math.log10(lfs_val)
    print(f"Free Space Path Loss (Lfs) at {freq_Hz/1e9} GHz over {dist_m/1e3} km = {Lfs_dB:.2f} dB")

    # Calculate Received Power at the input of the Rx antenna
    Prx_dBm = Ptx_dBm + Gtx_dB - Ltx_total_dB - Lfs_dB + Grx_dB
    print(f"\nReceived Power (Prx) = Ptx + Gtx - Ltx - Lfs + Grx")
    print(f"Prx = {Ptx_dBm} dBm + {Gtx_dB} dB - {Ltx_total_dB:.2f} dB - {Lfs_dB:.2f} dB + {Grx_dB} dB = {Prx_dBm:.2f} dBm\n")
    
    print("--- Step 2: Calculate System Noise Figure (NF_sys) and Input-Referred Noise (N_in) ---")

    # Calculate Cascaded Noise Figure using Friis Formula
    # F_sys = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
    F_sys = 0.0
    G_cascade_lin = 1.0
    print("Calculating system noise figure using Friis formula:")
    for i, stage in enumerate(rx_chain):
        F_stage = db_to_linear(stage['nf_db'])
        G_stage = db_to_linear(stage['gain_db'])
        
        contribution = (F_stage - 1) / G_cascade_lin
        if i == 0:
            contribution = F_stage # First stage contribution is just F1
        
        F_sys += contribution
        G_cascade_lin *= G_stage
        print(f"  - Stage {i+1} ({stage['name']}): F={F_stage:.2f}, Cumulative Gain={linear_to_db(G_cascade_lin):.1f} dB, Contribution to F_sys={contribution:.4f}")

    NF_sys_dB = linear_to_db(F_sys)
    print(f"\nTotal System Noise Factor (F_sys) = {F_sys:.4f}")
    print(f"Total System Noise Figure (NF_sys) = {NF_sys_dB:.2f} dB")

    # Calculate Input-Referred Noise Power (N_in)
    noise_floor_dbm_hz = 10 * math.log10(k * T * 1000) # -174 dBm/Hz
    noise_bw_factor_db = 10 * math.log10(B_Hz)
    
    N_in_dBm = noise_floor_dbm_hz + noise_bw_factor_db + NF_sys_dB
    
    print(f"\nInput-Referred Noise (N_in) = Thermal Noise Floor + 10*log10(Bandwidth) + NF_sys")
    print(f"N_in = {noise_floor_dbm_hz:.2f} dBm/Hz + {noise_bw_factor_db:.2f} dB-Hz + {NF_sys_dB:.2f} dB = {N_in_dBm:.2f} dBm\n")

    print("--- Step 3: Calculate Final Signal-to-Noise Ratio (SNR) ---")
    SNR_dB = Prx_dBm - N_in_dBm
    print(f"SNR (dB) = Prx (dBm) - N_in (dBm)")
    print(f"SNR = {Prx_dBm:.2f} dBm - ({N_in_dBm:.2f} dBm) = {SNR_dB:.2f} dB")
    
    return SNR_dB

if __name__ == '__main__':
    final_snr = main()
    # The final answer is wrapped in <<<>>>
    print(f"\n<<<Result>>>")
    print(f"<<<{final_snr:.2f}>>>")
