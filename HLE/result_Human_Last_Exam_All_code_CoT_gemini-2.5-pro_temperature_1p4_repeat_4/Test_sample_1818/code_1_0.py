import math

def calculate_toa(payload_bytes, sf, bw_khz=125, cr_code=1, preamble_len=8, header_enabled=True):
    """Calculates the Time on Air (ToA) for a LoRa packet."""
    bw_hz = bw_khz * 1000
    t_sym = (2**sf) / bw_hz  # Symbol time in seconds
    
    # Preamble duration
    t_preamble = (preamble_len + 4.25) * t_sym
    
    # Payload duration calculation
    h = 1 if header_enabled else 0
    de = 1 if sf >= 11 else 0 # Data rate optimization
    
    # Number of payload symbols
    # Formula from Semtech datasheets
    payload_num_bits = 8 * payload_bytes - 4 * sf + 28 + 16 - 20 * h
    numerator = max(0, payload_num_bits)
    denominator = 4 * (sf - 2 * de)
    
    n_payload_sym = 8 + math.ceil(numerator / denominator) * (cr_code + 4)
    
    t_payload = n_payload_sym * t_sym
    
    total_toa_s = t_preamble + t_payload
    return total_toa_s

def find_optimal_config(target_path_loss_db):
    """
    Finds the optimal SF and TP for a given target path loss.
    """
    # --- Parameters ---
    PAYLOAD_BYTES = 100
    BANDWIDTH_KHZ = 125
    CODING_RATE_CODE = 1  # 4/5
    
    SPREADING_FACTORS = list(range(7, 13)) # SF7 to SF12
    TRANSMIT_POWERS_DBM = list(range(2, 15, 2)) # 2dBm to 14dBm in 2dB steps
    
    # SNR required for demodulation at the gateway for each SF
    SNR_THRESHOLDS_DB = {
        7: -7.5, 8: -10, 9: -12.5, 10: -15, 11: -17.5, 12: -20
    }
    
    # For Rician K=3dB channel, a fading margin is needed for 99% reliability (1% PER)
    RICIAN_FADING_MARGIN_DB = 10.5 
    
    # Receiver noise floor calculation
    NOISE_FIGURE_DB = 6 # Typical gateway noise figure
    THERMAL_NOISE_DBM = -174 + 10 * math.log10(BANDWIDTH_KHZ * 1000)
    RECEIVER_NOISE_FLOOR_DBM = THERMAL_NOISE_DBM + NOISE_FIGURE_DB
    
    # --- Analysis ---
    best_config = None
    min_energy_uj = float('inf')

    print(f"Analyzing for a target path loss of {target_path_loss_db} dB...")
    print("-" * 50)
    print(f"{'SF':<4}{'Req TP (dBm)':<15}{'Actual TP (dBm)':<18}{'ToA (ms)':<12}{'Energy (uJ)':<12}")
    print("-" * 50)

    for sf in SPREADING_FACTORS:
        # Calculate total required SNR at the receiver
        required_avg_snr_db = SNR_THRESHOLDS_DB[sf] + RICIAN_FADING_MARGIN_DB
        
        # Calculate the required transmit power to overcome the path loss
        # PathLoss <= TP - NoiseFloor - RequiredSNR
        # TP >= PathLoss + NoiseFloor + RequiredSNR
        required_tp_dbm = target_path_loss_db + RECEIVER_NOISE_FLOOR_DBM + required_avg_snr_db

        # Find the lowest available TP that meets the requirement
        actual_tp_dbm = None
        for tp in TRANSMIT_POWERS_DBM:
            if tp >= required_tp_dbm:
                actual_tp_dbm = tp
                break
        
        # If no available TP can meet the requirement, this SF is not viable
        if actual_tp_dbm is None:
            print(f"{sf:<4}{required_tp_dbm:<15.1f}{'Not Possible':<18}{'-':<12}{'-':<12}")
            continue

        # Calculate Time on Air
        toa_s = calculate_toa(PAYLOAD_BYTES, sf, BANDWIDTH_KHZ, CODING_RATE_CODE)
        
        # Calculate energy consumption: E = P * t
        # Power (mW) = 10^(Power(dBm)/10)
        # Power (W) = 10^((Power(dBm)-30)/10)
        # Energy (uJ) = Power (mW) * ToA (ms)
        power_mw = 10**(actual_tp_dbm / 10)
        energy_uj = power_mw * (toa_s * 1000)

        print(f"{sf:<4}{required_tp_dbm:<15.1f}{actual_tp_dbm:<18}{toa_s * 1000:<12.2f}{energy_uj:<12.2f}")
        
        # Check if this configuration is the new best
        if energy_uj < min_energy_uj:
            min_energy_uj = energy_uj
            best_config = {
                "sf": sf,
                "tp": actual_tp_dbm,
                "toa_ms": toa_s * 1000,
                "energy_uj": energy_uj
            }
            
    # --- Final Result ---
    if best_config:
        print("\n" + "=" * 50)
        print("Optimal Configuration Found:")
        print(f"  Spreading Factor: SF{best_config['sf']}")
        print(f"  Transmission Power: {best_config['tp']} dBm")
        print(f"  Energy per Transmission: {best_config['energy_uj']:.2f} uJ")
        print("\nThis configuration minimizes energy for the assumed path loss while meeting the 1% PER requirement.")
        print("Final energy calculation:")
        
        power_mw = 10**(best_config['tp'] / 10)
        
        print(f"{best_config['energy_uj']:.2f} uJ = {best_config['toa_ms']:.2f} ms * {power_mw:.2f} mW (which is 10^({best_config['tp']}/10))")

    else:
        print("\n" + "=" * 50)
        print("No suitable configuration found for the given path loss.")
        
    return best_config

# Based on the analysis, a typical challenging urban path loss is assumed to find a single optimal point.
ASSUMED_PATH_LOSS_DB = 130
optimal_configuration = find_optimal_config(ASSUMED_PATH_LOSS_DB)
if optimal_configuration:
    final_sf = optimal_configuration["sf"]
    final_tp = optimal_configuration["tp"]

    # The final answer format is strict.
    # We output it here after all the explanatory text.
    print(f"\n<<<SF{final_sf}, {final_tp} dBm>>>")
