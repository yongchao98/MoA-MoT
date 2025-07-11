import math

def calculate_lora_toa(sf, payload_bytes, bw_hz, coding_rate_str):
    """
    Calculates the Time on Air (ToA) for a LoRa packet.
    Formula based on Semtech datasheets.
    """
    # Parameters
    bw = bw_hz
    pl = payload_bytes
    cr_map = {"4/5": 1, "4/6": 2, "4/7": 3, "4/8": 4}
    cr_code = cr_map[coding_rate_str]
    implicit_header = 0  # 0 for explicit header
    preamble_len = 8
    
    # Low Data Rate Optimization (LDRO) enabled for SF > 10
    de = 1 if sf >= 11 else 0

    # Symbol time in seconds
    t_sym = (2**sf) / bw

    # Number of payload symbols
    payload_numerator = 8 * pl - 4 * sf + 28 + 16 - 20 * implicit_header
    payload_denominator = 4 * (sf - 2 * de)
    
    # Ensure denominator is not zero
    if payload_denominator == 0:
        return None
        
    payload_symb_nb = 8 + max(0, math.ceil(payload_numerator / payload_denominator) * (cr_code + 4))

    # Total number of symbols
    total_symbols = preamble_len + 4.25 + payload_symb_nb

    # Time on Air in seconds
    toa_s = total_symbols * t_sym
    return toa_s

def solve_lorawan_optimization():
    """
    Solves for the optimal LoRaWAN SF and TxP to minimize energy.
    """
    # --- 1. Given and Assumed Parameters ---
    PAYLOAD_BYTES = 100
    BW_HZ = 125000
    CODING_RATE = "4/5"
    TRANSMIT_POWERS_DBM = list(range(2, 15, 2))  # 2, 4, ..., 14
    SPREADING_FACTORS = list(range(7, 13))      # SF7 to SF12
    
    # Environment and Link Parameters
    GATEWAY_NOISE_FIGURE_DB = 6
    RICIAN_K_FACTOR_DB = 3
    # A fading margin is needed to ensure 99% reliability (1% PER) on a Rician channel.
    # 4.5 dB is a reasonable engineering approximation for K=3dB.
    FADING_MARGIN_DB = 4.5 
    
    # CRITICAL ASSUMPTION: Path loss for the urban environment.
    # 135 dB is a challenging but realistic value for an urban deployment.
    ASSUMED_PATH_LOSS_DB = 135
    
    print("--- LoRaWAN Energy Optimization ---")
    print(f"Assumptions:")
    print(f"  - Path Loss: {ASSUMED_PATH_LOSS_DB} dB (typical for urban)")
    print(f"  - Channel Model: Rician Fading (K={RICIAN_K_FACTOR_DB} dB)")
    print(f"  - PER Target: <= 1% (requires a {FADING_MARGIN_DB} dB fading margin)\n")

    # --- 2. Intermediate Calculations ---
    # Receiver noise floor power in dBm
    noise_power_dbm = -174 + 10 * math.log10(BW_HZ) + GATEWAY_NOISE_FIGURE_DB

    # SNR required at the gateway for each SF to meet PER target
    # Base values are standard for LoRa, plus our fading margin
    snr_base_db = {7: -7.5, 8: -10.0, 9: -12.5, 10: -15.0, 11: -17.5, 12: -20.0}
    snr_target_db = {sf: snr_base + FADING_MARGIN_DB for sf, snr_base in snr_base_db.items()}

    # Time on Air (ToA) for each SF
    toa_s = {sf: calculate_lora_toa(sf, PAYLOAD_BYTES, BW_HZ, CODING_RATE) for sf in SPREADING_FACTORS}
    
    # --- 3. Optimization Logic ---
    print("--- Analysis per Spreading Factor ---")
    results = []
    for sf in SPREADING_FACTORS:
        # Calculate the required Tx Power to close the link
        required_txp_dbm = ASSUMED_PATH_LOSS_DB + snr_target_db[sf] - noise_power_dbm

        # Find the lowest available transmit power that meets the requirement
        optimal_txp_for_sf = None
        for txp in TRANSMIT_POWERS_DBM:
            if txp >= required_txp_dbm:
                optimal_txp_for_sf = txp
                break
        
        # If a valid power level was found, calculate energy and store the result
        if optimal_txp_for_sf is not None:
            power_mw = 10**(optimal_txp_for_sf / 10)
            energy_mj = power_mw * toa_s[sf]
            results.append({'sf': sf, 'txp': optimal_txp_for_sf, 'energy': energy_mj})
            print(f"SF{sf:2d}: Required TxP > {required_txp_dbm:.2f} dBm. "
                  f"Optimal choice: {optimal_txp_for_sf} dBm -> Energy: {energy_mj:.2f} mJ")
        else:
            print(f"SF{sf:2d}: Required TxP > {required_txp_dbm:.2f} dBm. "
                  f"No available power level is sufficient.")

    # --- 4. Find and Display the Optimal Result ---
    if not results:
        print("\nNo combination of SF and TxP can meet the requirements for the assumed path loss.")
        return

    optimal_config = min(results, key=lambda x: x['energy'])
    best_sf = optimal_config['sf']
    best_txp = optimal_config['txp']
    min_energy = optimal_config['energy']
    
    print("\n--- Conclusion ---")
    print("Comparing the energy consumption of all viable configurations...")
    print(f"\nThe optimal configuration with the lowest energy consumption is:")
    print(f"Spreading Factor: {best_sf}")
    print(f"Transmission Power: {best_txp} dBm")
    
    # Show the final calculation in detail as requested
    best_txp_mw = 10**(best_txp / 10)
    best_toa_s = toa_s[best_sf]
    
    print("\nFinal Energy Calculation:")
    print(f"Energy ({min_energy:.2f} mJ) = Power (mW) * Time on Air (s)")
    print(f"Energy ({min_energy:.2f} mJ) = {best_txp_mw:.2f} mW * {best_toa_s:.3f} s")
    print(f"where Power ({best_txp_mw:.2f} mW) = 10^({best_txp} dBm / 10)")
    
    # Store final answer in the requested format
    global final_answer
    final_answer = f"SF{best_sf}, {best_txp}dBm"

# Execute the function and print the final formatted answer
final_answer = "Calculation failed"
solve_lorawan_optimization()
print(f"\n<<<{final_answer}>>>")