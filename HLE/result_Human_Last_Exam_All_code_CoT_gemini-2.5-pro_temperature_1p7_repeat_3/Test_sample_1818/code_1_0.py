import math

def calculate_lora_toa(sf, payload_bytes, bw=125000, n_preamble=8, explicit_header=True, crc=True, coding_rate_den=5):
    """
    Calculates the Time on Air for a LoRa packet.
    """
    # Symbol time
    t_symbol = (2**sf) / bw

    # Preamble duration
    t_preamble = (n_preamble + 4.25) * t_symbol

    # Low Data Rate Optimization (for SF11 and SF12 on 125kHz)
    de = 0
    if sf >= 11:
        de = 1
    
    # Header configuration
    h = 0 if explicit_header else 1
    
    # CRC configuration
    crc_val = 1 if crc else 0
    
    # Coding rate value for formula
    cr = max(coding_rate_den - 4, 1)

    # Number of payload symbols
    payload_symbols_numerator = 8 * payload_bytes - 4 * sf + 28 + 16 * crc_val - 20 * h
    payload_symbols_denominator = 4 * (sf - 2 * de)
    
    if payload_symbols_denominator == 0:
        return None, None # Avoid division by zero

    payload_symbols = 8 + max(0, math.ceil(payload_symbols_numerator / payload_symbols_denominator) * cr)

    # Payload duration
    t_payload = payload_symbols * t_symbol
    
    # Total ToA
    total_toa_s = t_preamble + t_payload
    return total_toa_s * 1000  # Return in milliseconds

def solve_lorawan_energy():
    """
    Solves for the optimal SF and TxP to minimize energy consumption.
    """
    # --- Given Parameters ---
    PAYLOAD_BYTES = 100
    BANDWIDTH = 125000  # 125 kHz
    CODING_RATE = 4/5
    AVAILABLE_TX_POWER_DBM = list(range(2, 15, 2))
    SPREADING_FACTORS = list(range(7, 13))

    # --- Modeling & Assumptions ---
    # SNR thresholds for demodulation at near-zero PER
    SNR_THRESHOLDS = {7: -7.5, 8: -10.0, 9: -12.5, 10: -15.0, 11: -17.5, 12: -20.0}
    # Fading margin for Rician K=3dB at 99% reliability (1% PER)
    FADING_MARGIN_DB = 5.0
    # Assumed scenario for a concrete answer
    ASSUMED_PATH_LOSS_DB = 130.0
    ASSUMED_GATEWAY_NOISE_FLOOR_DBM = -117.0
    
    print("LoRaWAN Energy Optimization Analysis\n")
    print("--- Assumptions ---")
    print(f"Payload: {PAYLOAD_BYTES} bytes, Bandwidth: {BANDWIDTH/1000} kHz, Coding Rate: {int(CODING_RATE*5)}/5")
    print(f"Rician Fading Margin (K=3dB, PER<=1%): {FADING_MARGIN_DB} dB")
    print(f"Assumed Urban Path Loss: {ASSUMED_PATH_LOSS_DB} dB")
    print(f"Assumed Gateway Noise Floor: {ASSUMED_GATEWAY_NOISE_FLOOR_DBM} dBm\n")

    results = []
    
    print("--- Calculations ---")
    print("-" * 105)
    print(f"{'SF':<5} | {'ToA (ms)':<12} | {'SNR_req (dB)':<14} | {'TxP_req (dBm)':<15} | {'TxP_set (dBm)':<15} | {'Power (mW)':<12} | {'Energy (mJ)':<15}")
    print("-" * 105)

    for sf in SPREADING_FACTORS:
        # 1. Calculate Time on Air
        toa_ms = calculate_lora_toa(sf, PAYLOAD_BYTES, bw=BANDWIDTH, coding_rate_den=5)
        
        # 2. Determine Required SNR
        snr_req_db = SNR_THRESHOLDS[sf] + FADING_MARGIN_DB
        
        # 3. Determine Required Transmit Power for the assumed scenario
        # TxP_req = SNR_req + PathLoss - (G_tx + G_rx) + NoiseFloor
        # We simplify G_tx and G_rx to 0. So: TxP_req = SNR_req + PathLoss + NoiseFloor
        txp_req_dbm = snr_req_db + ASSUMED_PATH_LOSS_DB + ASSUMED_GATEWAY_NOISE_FLOOR_DBM

        # 4. Find the minimum available Tx Power that meets the requirement
        txp_assigned_dbm = None
        for p in AVAILABLE_TX_POWER_DBM:
            if p >= txp_req_dbm:
                txp_assigned_dbm = p
                break
        
        # 5. Calculate energy if a valid power level was found
        if txp_assigned_dbm is not None:
            power_mw = 10**(txp_assigned_dbm / 10)
            energy_mj = toa_ms * power_mw / 1000 # Convert mW*ms to mJ
            results.append({
                'sf': sf, 'toa_ms': toa_ms, 'snr_req_db': snr_req_db,
                'txp_req_dbm': txp_req_dbm, 'txp_assigned_dbm': txp_assigned_dbm,
                'power_mw': power_mw, 'energy_mj': energy_mj
            })
            print(f"SF{sf:<2} | {toa_ms:<12.1f} | {snr_req_db:<14.1f} | {txp_req_dbm:<15.1f} | {txp_assigned_dbm:<15} | {power_mw:<12.2f} | {energy_mj:<15.1f}")
        else:
            print(f"SF{sf:<2} | {toa_ms:<12.1f} | {snr_req_db:<14.1f} | {txp_req_dbm:<15.1f} | {'(>14 dBm)':<15} | {'N/A':<12} | {'N/A':<15}")

    print("-" * 105)

    # 6. Find the optimal configuration with minimum energy
    if not results:
        print("\nNo configuration could meet the requirements with the available power levels.")
        return

    optimal_config = min(results, key=lambda x: x['energy_mj'])

    print("\n--- Conclusion ---")
    print("To minimize energy consumption for the given scenario, the network server should select the configuration with the lowest energy cost.")
    
    sf_opt = optimal_config['sf']
    txp_opt_dbm = optimal_config['txp_assigned_dbm']
    toa_opt = optimal_config['toa_ms']
    power_opt_mw = optimal_config['power_mw']
    energy_opt_mj = optimal_config['energy_mj']

    print(f"\nThe optimal configuration is SF{sf_opt} with a transmit power of {txp_opt_dbm} dBm.")
    print("This results in the lowest energy consumption while maintaining a Packet Error Rate below 1%.")
    print("\nFinal Optimal Energy Equation:")
    print(f"Energy ({energy_opt_mj:.1f} mJ) = Time on Air ({toa_opt:.1f} ms) * Transmit Power ({txp_opt_dbm} dBm or {power_opt_mw:.2f} mW)")
    
    # Final answer in the requested format
    print(f"\n<<<{sf_opt}, {txp_opt_dbm} dBm>>>")

solve_lorawan_energy()