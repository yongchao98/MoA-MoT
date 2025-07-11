import math

def calculate_optimal_settings():
    """
    Calculates the optimal LoRaWAN Spreading Factor and Transmit Power
    to minimize energy consumption for the given scenario.
    """

    # --- Given Parameters ---
    PAYLOAD_BYTES = 100
    BANDWIDTH_HZ = 125000
    CODING_RATE_CODE = 1  # For CR = 4/5
    PREAMBLE_SYMBOLS = 8
    HEADER_ENABLED = 0  # 0 for explicit header
    CRC_ENABLED = 1  # 1 for CRC on
    AVAILABLE_TX_POWERS_DBM = list(range(2, 15, 2))

    # SNR required for ~1% PER for each SF
    SNR_REQUIREMENTS_DB = {
        7: -7.5,
        8: -10,
        9: -12.5,
        10: -15,
        11: -17.5,
        12: -20
    }

    # Assumed constant for link loss and noise floor based on "urban environment"
    # L_const = Path Loss + Noise Figure - Antenna Gains
    # A value of 20 dB implies a significant path loss (~140-145 dB).
    L_CONST_DB = 20

    print(f"Analyzing energy consumption for a {PAYLOAD_BYTES}-byte payload...")
    print(f"Assuming a constant link loss (L_const) of {L_CONST_DB} dB.\n")

    results = []

    # Iterate through each Spreading Factor
    for sf in range(7, 13):
        # 1. Calculate Time on Air (ToA)
        is_low_dr = 1 if sf >= 11 else 0
        symbol_time_s = (2**sf) / BANDWIDTH_HZ

        payload_num_part = 8 * PAYLOAD_BYTES - 4 * sf + 28 + 16 * CRC_ENABLED - 20 * HEADER_ENABLED
        payload_den_part = 4 * (sf - 2 * is_low_dr)
        
        # Ceil division and multiplication by coding rate factor
        num_payload_symbols = 8 + max(0, math.ceil(payload_num_part / payload_den_part) * (CODING_RATE_CODE + 4))

        toa_s = (PREAMBLE_SYMBOLS + 4.25 + num_payload_symbols) * symbol_time_s

        # 2. Determine Required Transmit Power
        snr_req = SNR_REQUIREMENTS_DB[sf]
        tx_p_req = snr_req + L_CONST_DB

        # Find the smallest available power that meets the requirement
        tx_p_actual = None
        for p in AVAILABLE_TX_POWERS_DBM:
            if p >= tx_p_req:
                tx_p_actual = p
                break
        
        print(f"--- SF{sf} ---")
        print(f"   - Time on Air (ToA): {toa_s * 1000:.2f} ms")
        print(f"   - Required SNR: {snr_req} dB")
        print(f"   - Required Tx Power: {snr_req} dB (SNR) + {L_CONST_DB} dB (Loss) = {tx_p_req:.2f} dBm")

        if tx_p_actual is None:
            print("   - No available power level is sufficient for this SF.\n")
            continue
            
        print(f"   - Assigned Tx Power: {tx_p_actual} dBm")
        
        # 3. Calculate Relative Energy
        tx_p_linear = 10**(tx_p_actual / 10)
        energy = tx_p_linear * toa_s
        
        results.append({'sf': sf, 'tx_p': tx_p_actual, 'energy': energy})

        print(f"   - Relative Energy Calculation:")
        print(f"     10^({tx_p_actual} / 10) * {toa_s:.4f} = {energy:.4f}")
        print("")


    # 4. Find the optimal setting with the minimum energy
    if not results:
        print("Could not find a valid setting for the given path loss.")
        return

    optimal_result = min(results, key=lambda x: x['energy'])
    
    print("--- Conclusion ---")
    print("The optimal setting is the one with the lowest relative energy consumption.")
    print(f"Optimal Spreading Factor: SF{optimal_result['sf']}")
    print(f"Optimal Transmission Power: {optimal_result['tx_p']} dBm")


calculate_optimal_settings()
<<<SF8, 10 dBm>>>