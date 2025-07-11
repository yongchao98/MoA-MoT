import numpy as np
import math

def calculate_toa(payload_len, sf, bw=125000, cr_code=1, n_preamble=8, crc=True, header=True):
    """
    Calculates the Time on Air (ToA) for a LoRa packet.
    """
    # Low Data Rate Optimization is used for SF11 and SF12 with 125kHz BW
    de = 1 if sf >= 11 else 0
    
    # Symbol time
    t_sym = (2**sf) / bw
    
    # Preamble duration
    t_preamble = (n_preamble + 4.25) * t_sym
    
    # Number of payload symbols
    # Formula from Semtech datasheets
    # 8*PL - 4*SF + 28 + 16(CRC) - 20*H
    # H=0 for header enabled, H=1 for headerless
    h = 0 if header else 1
    crc_val = 16 if crc else 0
    
    numerator = 8 * payload_len - 4 * sf + 28 + crc_val - 20 * h
    denominator = 4 * (sf - 2 * de)
    
    n_payload_sym = 8 + max(math.ceil(numerator / denominator) * (cr_code + 4), 0)
    
    # Payload duration
    t_payload = n_payload_sym * t_sym
    
    # Total ToA in seconds
    toa_s = t_preamble + t_payload
    return toa_s * 1000 # Return in milliseconds

def solve_lora_optimization():
    """
    Finds the optimal SF and TxP for the given LoRaWAN scenario.
    """
    # --- Given Parameters ---
    payload_len = 100
    coding_rate_denom = 5 # CR = 4/5, so CR_code = 1
    cr_code = coding_rate_denom - 4
    bandwidth = 125000
    tx_powers_dbm = list(range(2, 15, 2))
    spreading_factors = list(range(7, 13))
    
    # --- Assumptions ---
    # Assume a path loss typical for a challenging urban environment (cell edge)
    path_loss_db = 148.0
    # Typical noise floor for a 125kHz channel
    noise_floor_dbm = -126.0
    
    # --- Fading Margin Calculation ---
    # For a Rician channel with K=3dB, a 99% success rate (1% PER) requires
    # a fading margin of approximately 6.54 dB.
    fading_margin_db = 6.54
    
    # SNR requirements (demodulation floor) for each SF
    snr_base_req = {
        7: -7.5, 8: -10.0, 9: -12.5,
        10: -15.0, 11: -17.5, 12: -20.0
    }
    
    print("--- LoRaWAN Energy Optimization ---")
    print(f"Assumed Path Loss: {path_loss_db} dB")
    print(f"Assumed Noise Floor: {noise_floor_dbm} dBm")
    print(f"Required Fading Margin for PER < 1% (Rician K=3dB): {fading_margin_db} dB\n")
    
    results = []
    
    for sf in spreading_factors:
        # 1. Calculate Time on Air (ToA)
        toa_ms = calculate_toa(payload_len, sf, cr_code=cr_code)
        
        # 2. Calculate required received power at the gateway
        total_snr_req = snr_base_req[sf] + fading_margin_db
        req_rx_power_dbm = noise_floor_dbm + total_snr_req
        
        # 3. Calculate required transmit power from the device
        req_tx_power_dbm = path_loss_db + req_rx_power_dbm
        
        # 4. Find the smallest available TxP that meets the requirement
        actual_tx_power_dbm = -1
        for txp in tx_powers_dbm:
            if txp >= req_tx_power_dbm:
                actual_tx_power_dbm = txp
                break
        
        print(f"--- Analyzing SF{sf} ---")
        print(f"  Time on Air (ToA): {toa_ms:.2f} ms")
        print(f"  Required Rx Power: {req_rx_power_dbm:.2f} dBm")
        print(f"  Required Tx Power: {req_tx_power_dbm:.2f} dBm")
        
        # 5. If a viable TxP is found, calculate the energy index
        if actual_tx_power_dbm != -1:
            tx_power_mw = 10**(actual_tx_power_dbm / 10)
            energy_index = toa_ms * tx_power_mw
            results.append({
                'sf': sf,
                'txp_dbm': actual_tx_power_dbm,
                'energy': energy_index,
                'toa_ms': toa_ms,
                'txp_mw': tx_power_mw
            })
            print(f"  Assigned Tx Power: {actual_tx_power_dbm} dBm (Viable)")
            print(f"  Relative Energy Index: {energy_index:.2f}")
        else:
            print(f"  Assigned Tx Power: Not possible (Max power is {max(tx_powers_dbm)} dBm)")
        print("-" * 25)

    # 6. Find the configuration with the minimum energy index
    if not results:
        print("\nNo viable SF/TxP combination found for the assumed path loss.")
        return
        
    optimal_config = min(results, key=lambda x: x['energy'])
    
    print("\n--- Optimal Configuration Found ---")
    print(f"Optimal Spreading Factor: SF{optimal_config['sf']}")
    print(f"Optimal Transmission Power: {optimal_config['txp_dbm']} dBm")
    
    print("\nThis configuration provides the required link reliability with the lowest energy consumption.")
    print("Final energy calculation for the optimal setting:")
    
    # Output each number in the final equation as requested
    toa = optimal_config['toa_ms']
    txp_dbm = optimal_config['txp_dbm']
    txp_mw = optimal_config['txp_mw']
    energy = optimal_config['energy']
    
    print(f"Energy Index = ToA * 10^(TxP / 10)")
    print(f"             = {toa:.2f} * 10^({txp_dbm} / 10)")
    print(f"             = {toa:.2f} * {txp_mw:.2f}")
    print(f"             = {energy:.2f}")

    print(f"\n<<<SF{optimal_config['sf']} at {optimal_config['txp_dbm']} dBm>>>")

solve_lora_optimization()