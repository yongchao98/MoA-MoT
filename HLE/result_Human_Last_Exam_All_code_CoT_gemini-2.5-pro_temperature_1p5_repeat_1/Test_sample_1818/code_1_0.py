import math

def calculate_lora_toa(sf, payload_bytes, bw=125000, cr_code=1, n_preamble=8, header_enabled=True, crc_enabled=True):
    """
    Calculates the Time on Air (ToA) for a LoRaWAN packet.
    """
    # Low Data Rate Optimization is enabled for SF11 and SF12 in EU868
    low_dr_opt = 1 if sf >= 11 else 0
    
    # Symbol time
    t_sym = (2**sf) / bw
    
    # Preamble duration
    t_preamble = (n_preamble + 4.25) * t_sym
    
    # Number of payload symbols
    payload_symbol_numerator = (8 * payload_bytes - 4 * sf + 28 + (16 if crc_enabled else 0) - (20 if header_enabled else 0))
    payload_symbol_denominator = 4 * (sf - 2 * low_dr_opt)
    
    n_payload_sym = 8 + max(0, math.ceil(payload_symbol_numerator / payload_symbol_denominator) * (cr_code + 4))
    
    # Payload duration
    t_payload = n_payload_sym * t_sym
    
    # Total ToA
    toa = t_preamble + t_payload
    return toa

def solve_lorawan_optimization():
    """
    Finds the optimal SF and TP for the given LoRaWAN scenario.
    """
    # --- Step 1: Define Parameters ---
    payload_size = 100  # bytes
    lorawan_header_size = 13 # bytes
    total_payload = payload_size + lorawan_header_size

    spreading_factors = list(range(7, 13))  # SF7 to SF12
    transmit_powers_dbm = list(range(2, 15, 2))  # 2 to 14 dBm in 2 dB steps
    
    # SNR thresholds for demodulation in AWGN channel
    snr_thresholds_awgn = {7: -7.5, 8: -10, 9: -12.5, 10: -15, 11: -17.5, 12: -20}
    
    # Fading margin for Rician channel (K=3dB) to achieve PER <= 1%
    fading_margin_db = 4.5
    
    # Assumed target Link Gain for an urban environment
    target_link_gain_db = 22.0

    # --- Step 2: Calculate performance for all settings ---
    configs = []
    for sf in spreading_factors:
        # Calculate ToA (depends on SF only)
        toa_s = calculate_lora_toa(sf=sf, payload_bytes=total_payload, cr_code=1) # CR 4/5 -> code 1
        
        # Calculate required SNR including fading margin
        snr_req_db = snr_thresholds_awgn[sf] + fading_margin_db
        
        for tp_dbm in transmit_powers_dbm:
            # Convert power to Watts for energy calculation
            tp_watts = 10**((tp_dbm - 30) / 10)
            
            # Calculate energy in Joules
            energy_j = tp_watts * toa_s
            
            # Calculate the Link Gain this configuration supports
            link_gain_db = tp_dbm - snr_req_db
            
            configs.append({
                "sf": sf,
                "tp_dbm": tp_dbm,
                "toa_s": toa_s,
                "energy_mj": energy_j * 1000,
                "link_gain_db": link_gain_db,
                "snr_req_db": snr_req_db,
                "tp_watts": tp_watts
            })

    # --- Step 3: Find the efficiency frontier (Pareto optimal points) ---
    frontier = []
    for candidate in sorted(configs, key=lambda x: x['link_gain_db']):
        is_dominated = False
        for other in configs:
            if other == candidate:
                continue
            # A candidate is dominated if another config offers more (or equal) gain
            # for less (or equal) energy.
            if other['link_gain_db'] >= candidate['link_gain_db'] and other['energy_mj'] <= candidate['energy_mj']:
                 # Strict inequality check to avoid self-domination on ties
                if other['link_gain_db'] > candidate['link_gain_db'] or other['energy_mj'] < candidate['energy_mj']:
                    is_dominated = True
                    break
        if not is_dominated:
            frontier.append(candidate)

    # --- Step 4: Select the best option from the frontier for the urban scenario ---
    best_option = None
    min_energy = float('inf')

    for config in sorted(frontier, key=lambda x: x['energy_mj']):
        if config['link_gain_db'] >= target_link_gain_db:
            if config['energy_mj'] < min_energy:
                min_energy = config['energy_mj']
                best_option = config

    # --- Step 5: Print the results and the final equation ---
    if best_option:
        print(f"Optimal Configuration for Urban Environment (Target Link Gain >= {target_link_gain_db} dB):")
        print("-" * 60)
        print(f"Spreading Factor (SF): {best_option['sf']}")
        print(f"Transmission Power (TP): {best_option['tp_dbm']} dBm")
        print("\nThis configuration results in:")
        print(f"  - Time on Air: {best_option['toa_s']:.3f} s")
        print(f"  - Energy per Transmission: {best_option['energy_mj']:.2f} mJ")
        print(f"  - Achieved Link Gain: {best_option['link_gain_db']:.1f} dB")
        
        print("\nFinal Calculation:")
        
        # Energy Calculation Output
        print(f"  Energy ({best_option['energy_mj']:.2f} mJ) = Power ({best_option['tp_watts']*1000:.2f} mW) * Time on Air ({best_option['toa_s']:.3f} s)")

        # Link Gain Calculation Output
        print(f"  Link Gain ({best_option['link_gain_db']:.1f} dB) = Transmit Power ({best_option['tp_dbm']} dBm) - Required SNR ({best_option['snr_req_db']:.1f} dB)")

    else:
        print(f"No available configuration meets the target link gain of {target_link_gain_db} dB.")

# Run the solver
solve_lorawan_optimization()