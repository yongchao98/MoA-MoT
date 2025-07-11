import math

def solve_lorawan_energy():
    """
    Calculates the optimal Spreading Factor and Transmit Power for a LoRaWAN device
    to minimize energy while maintaining a target PER.
    """
    # 1. Define constants from the problem description
    PAYLOAD_BYTES = 100
    BW = 125_000  # Bandwidth in Hz
    CR_CODE = 1   # Coding Rate 4/5 -> code = 1
    TRANSMIT_POWERS_DBM = list(range(2, 15, 2))  # 2 to 14 dBm in 2dB steps
    SPREADING_FACTORS = list(range(7, 13))  # SF7 to SF12

    # --- Assumed and Calculated Link Parameters ---
    # K-factor for Rician fading
    K_FACTOR_DB = 3
    # Target Packet Error Rate
    PER_TARGET = 0.01
    # For a Rician channel with K=3dB, a fading margin is needed to ensure PER < 1%.
    # This value is based on standard radio channel modeling.
    FADING_MARGIN_DB = 4.5
    
    # Gateway characteristics
    NOISE_FIGURE_DB = 6
    THERMAL_NOISE_DBM_PER_HZ = -174
    
    # A representative path loss for an urban environment
    PATH_LOSS_DB = 132

    # Standard LoRa preamble length
    PREAMBLE_SYMBOLS = 8
    
    # Calculate total noise power over the bandwidth
    noise_power_dbm = THERMAL_NOISE_DBM_PER_HZ + 10 * math.log10(BW) + NOISE_FIGURE_DB
    
    # SNR required for demodulation at each SF (standard values)
    SNR_DEMOD_THRESHOLDS_DB = {
        7: -7.5,
        8: -10,
        9: -12.5,
        10: -15,
        11: -17.5,
        12: -20,
    }

    # Helper function to calculate Time on Air (ToA)
    def calculate_toa_s(sf, payload_bytes):
        # Symbol Time in seconds
        ts = (2**sf) / BW
        
        # Preamble Time
        t_preamble = (PREAMBLE_SYMBOLS + 4.25) * ts
        
        # Payload Time calculation
        # Low Data Rate Optimization is enabled for SF11 and SF12 in EU868
        de = 1 if sf >= 11 else 0
        # Header is explicit, CRC is on
        ih = 0
        crc_on = 16

        payload_numerator = 8 * payload_bytes - 4 * sf + 28 + crc_on - 20 * ih
        payload_denominator = 4 * (sf - 2 * de)
        
        # Number of symbols in the payload
        symbols_payload = 8 + max(0, math.ceil(payload_numerator / payload_denominator) * (CR_CODE + 4))
        
        t_payload = symbols_payload * ts
        
        return t_preamble + t_payload

    # Store results
    results = []
    
    print(f"Analyzing for a representative path loss of {PATH_LOSS_DB} dB...")
    print("-" * 65)
    print(f"{'SF':<5} {'Req. Sens.':<12} {'Req. TP':<10} {'Asgn. TP':<10} {'ToA (s)':<10} {'Rel. Energy':<12}")
    print("-" * 65)

    # 2. Iterate through each SF to find the best configuration
    for sf in SPREADING_FACTORS:
        # a. Determine required receiver sensitivity
        snr_demod = SNR_DEMOD_THRESHOLDS_DB[sf]
        # Required average SNR includes the fading margin for 1% PER
        required_avg_snr = snr_demod + FADING_MARGIN_DB
        required_sensitivity_dbm = noise_power_dbm + required_avg_snr

        # b. Calculate minimum TP required to overcome path loss
        required_tp_dbm = required_sensitivity_dbm + PATH_LOSS_DB

        # c. Find the lowest available TP that meets the requirement
        assigned_tp_dbm = None
        for tp in TRANSMIT_POWERS_DBM:
            if tp >= required_tp_dbm:
                assigned_tp_dbm = tp
                break
        
        # If no available TP can close the link, this SF is not viable
        if assigned_tp_dbm is None:
            print(f"{sf:<5} {required_sensitivity_dbm:<-12.2f} {required_tp_dbm:<10.2f} {'N/A':<10} {'N/A':<10} {'Not Viable':<12}")
            continue

        # d. Calculate Time on Air for the payload
        toa_s = calculate_toa_s(sf, PAYLOAD_BYTES)

        # e. Compute relative energy (Power in mW * time in s)
        transmit_power_mw = 10**(assigned_tp_dbm / 10)
        relative_energy = transmit_power_mw * toa_s

        results.append({
            'sf': sf,
            'assigned_tp_dbm': assigned_tp_dbm,
            'relative_energy': relative_energy,
            'toa_s': toa_s
        })
        
        # Print results for this SF
        print(f"{sf:<5} {required_sensitivity_dbm:<12.2f} {required_tp_dbm:<10.2f} {assigned_tp_dbm:<10.1f} {toa_s:<10.3f} {relative_energy:<12.2f}")

    # 3. Identify the setting with the minimum energy consumption
    if not results:
        print("\nNo viable SF/TP combination found for the given path loss.")
        return

    optimal_setting = min(results, key=lambda x: x['relative_energy'])
    
    print("-" * 65)
    print("\nConclusion:")
    print("The optimal setting is the one that provides the required reliability (PER < 1%) with the lowest energy consumption.")
    print(f"Based on the analysis, the most energy-efficient configuration is:")
    
    opt_sf = optimal_setting['sf']
    opt_tp = optimal_setting['assigned_tp_dbm']
    opt_toa = optimal_setting['toa_s']
    opt_energy = optimal_setting['relative_energy']
    
    print(f"\nSpreading Factor: SF{opt_sf}")
    print(f"Transmission Power: {opt_tp} dBm")
    
    # Final equation as requested
    print(f"\nThis choice is based on the final energy calculation:")
    print(f"Energy = 10^({opt_tp} / 10) * {opt_toa:.4f} = {opt_energy:.3f} (relative units: mW * s)")

    return f"SF{opt_sf}, {opt_tp} dBm"


# Execute the solver and capture the final answer string
final_answer = solve_lorawan_energy()
print(f"\n<<<{final_answer}>>>")
