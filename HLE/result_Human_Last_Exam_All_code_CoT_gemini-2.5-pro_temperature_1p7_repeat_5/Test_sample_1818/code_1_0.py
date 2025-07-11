import math

def solve_lorawan_optimization():
    """
    Calculates the optimal Spreading Factor and Transmit Power for a LoRaWAN device
    to minimize energy consumption while meeting a PER target.
    """
    # --- 1. Define given and assumed parameters ---

    # LoRaWAN parameters
    payload_len = 100  # bytes
    bw = 125000  # Hz
    coding_rate_val = 4/5
    preamble_symbols = 8
    has_crc = True
    explicit_header = True
    available_tx_power_dbm = list(range(2, 15, 2))  # 2 to 14 dBm in 2 dB steps

    # Channel and Link Budget parameters
    path_loss_db = 135  # Assumed representative urban path loss
    k_factor_db = 3
    fading_margin_db = 3  # Fading margin for Rician K=3dB at 99% reliability
    gateway_noise_figure_db = 4  # Typical gateway noise figure
    thermal_noise_dbm = -174 + 10 * math.log10(bw)
    noise_floor_dbm = thermal_noise_dbm + gateway_noise_figure_db

    # SNR required for 1% PER in an AWGN channel (baseline)
    snr_req_awgn_db = {
        7: -6.0,
        8: -9.0,
        9: -12.0,
        10: -15.0,
        11: -17.5,
        12: -20.0,
    }

    results = []
    
    print("--- Analysis Report ---")
    print(f"Assumed Path Loss: {path_loss_db} dB")
    print(f"Required PER: <= 1%")
    print(f"Fading Margin for Rician (K={k_factor_db}dB) Channel: {fading_margin_db} dB\n")

    # --- 2. Iterate through each Spreading Factor ---
    for sf in range(7, 13):
        # --- a. Calculate Time on Air (ToA) ---
        de = 1 if sf >= 11 else 0  # Low Data Rate Optimization
        h = 0 if explicit_header else 1
        crc = 1 if has_crc else 0
        cr_code = 1 # for 4/5
        
        # Calculate number of payload symbols
        # Using the standard formula for LoRa packet structure
        payload_part_numerator = 8 * payload_len - 4 * sf + 28 + 16 * crc - 20 * h
        payload_part_denominator = 4 * (sf - 2 * de)
        
        num_payload_symbols = 8 + max(0, math.ceil(payload_part_numerator / payload_part_denominator) * (cr_code + 4))

        symbol_duration_s = (2**sf) / bw
        preamble_duration_s = (preamble_symbols + 4.25) * symbol_duration_s
        payload_duration_s = num_payload_symbols * symbol_duration_s
        toa_s = preamble_duration_s + payload_duration_s

        # --- b. Calculate Required Transmit Power ---
        # Add fading margin to AWGN SNR requirement
        snr_req_fading_db = snr_req_awgn_db[sf] + fading_margin_db
        
        # Calculate minimum required Tx Power using the link budget equation
        # Required_SNR = TP - PathLoss - NoiseFloor
        required_tp_dbm = snr_req_fading_db + path_loss_db + noise_floor_dbm

        # Find the lowest available TP that satisfies the requirement
        chosen_tp_dbm = -1
        for power in available_tx_power_dbm:
            if power >= required_tp_dbm:
                chosen_tp_dbm = power
                break
        
        # Skip if no available power level is sufficient
        if chosen_tp_dbm == -1:
            print(f"SF{sf:<3} | Infeasible for this path loss (Requires >{required_tp_dbm:.1f} dBm)")
            continue
            
        # --- c. Calculate Energy Consumption ---
        tx_power_watts = 0.001 * (10**(chosen_tp_dbm / 10))
        energy_joules = tx_power_watts * toa_s

        results.append({
            'sf': sf,
            'toa_ms': toa_s * 1000,
            'required_tp_dbm': required_tp_dbm,
            'chosen_tp_dbm': chosen_tp_dbm,
            'energy_mJ': energy_joules * 1000,
            'tx_power_watts': tx_power_watts
        })
        
        print(f"SF{sf:<3} | ToA: {toa_s*1000:7.1f} ms | Required TP: {required_tp_dbm:5.1f} dBm -> Chosen TP: {chosen_tp_dbm:2} dBm | Energy: {energy_joules*1000:5.2f} mJ")


    # --- 3. Find the optimal configuration (minimum energy) ---
    if not results:
        print("\nNo viable configuration found for the assumed path loss.")
        return

    optimal_config = min(results, key=lambda x: x['energy_mJ'])
    
    # --- 4. Print the final answer ---
    print("\n--- Optimal Configuration ---")
    print(f"The configuration with the lowest energy consumption is:")
    print(f"Spreading Factor: SF{optimal_config['sf']}")
    print(f"Transmission Power: {optimal_config['chosen_tp_dbm']} dBm")
    
    print("\nThis configuration results from the following calculation:")
    final_sf = optimal_config['sf']
    final_tp = optimal_config['chosen_tp_dbm']
    final_toa_s = optimal_config['toa_ms'] / 1000
    final_power_w = optimal_config['tx_power_watts']
    final_energy_j = optimal_config['energy_mJ'] / 1000
    
    print(f"   Time on Air for SF{final_sf}: {final_toa_s:.4f} s")
    print(f"   Transmit Power for {final_tp} dBm: 10^({final_tp}/10) * 0.001 = {final_power_w:.6f} W")
    print("\n   Final Energy Equation:")
    print(f"   Energy = Power * Time on Air")
    print(f"   Energy = {final_power_w:.6f} W * {final_toa_s:.4f} s = {final_energy_j:.6f} Joules")
    print(f"   Optimal Energy: {optimal_config['energy_mJ']:.3f} mJ")
    
    
solve_lorawan_optimization()
<<<SF10 at 4 dBm>>>