import math

def solve_lorawan_optimization():
    """
    This function determines the optimal LoRaWAN Spreading Factor and Transmission Power
    to minimize energy consumption for a given scenario.
    """

    # --- Plan Explanation ---
    # The goal is to find the optimal Spreading Factor (SF) and Transmission Power (TP)
    # for a LoRaWAN device in an urban environment to minimize energy while ensuring a
    # Packet Error Rate (PER) below 1%. The problem is under-specified without a
    # fixed path loss value. To find a single optimal solution, we will assume a
    # challenging but realistic path loss of 140 dB, typical for urban NLOS (Non-Line-of-Sight)
    # conditions. The optimal setting will be the one that can overcome this path loss
    # with the absolute minimum energy consumption.
    #
    # The steps are:
    # 1.  Define LoRaWAN and environmental parameters. We use a 3 dB link margin over the
    #     base SNR requirement to account for the Rician fading channel and meet the 1% PER.
    # 2.  Calculate the receiver sensitivity for each SF.
    # 3.  For the assumed 140 dB path loss, iterate through all (SF, TP) combinations to find
    #     which ones can successfully establish a link.
    # 4.  For each successful combination, calculate the total energy consumed per transmission.
    # 5.  Select the combination with the lowest calculated energy as the optimal setting.
    # 6.  Present the results, including the final energy calculation with all numerical values.

    # Step 1: Define parameters
    payload_bytes = 100
    bw = 125000  # Hz
    coding_rate_val = 1  # 1 corresponds to a Coding Rate of 4/5
    header_enabled = True
    crc_enabled = True
    preamble_symbols = 8
    
    sfs = list(range(7, 13))
    tps_dbm_available = list(range(2, 15, 2))

    assumed_path_loss_db = 140.0
    fading_margin_db = 3.0  # Margin for Rician K=3dB to achieve PER < 1%
    receiver_noise_figure_db = 6.0
    
    # Step 2: Calculate receiver sensitivity for each SF
    base_snr_req = {7: -7.5, 8: -10, 9: -12.5, 10: -15, 11: -17.5, 12: -20}
    noise_floor_dbm = -174 + 10 * math.log10(bw) + receiver_noise_figure_db
    
    # The required SNR is the base demodulation SNR plus the fading margin
    snr_req_with_margin = {sf: base_snr_req[sf] + fading_margin_db for sf in sfs}
    sensitivity = {sf: noise_floor_dbm + snr_req_with_margin[sf] for sf in sfs}

    valid_configs = []
    
    for sf in sfs:
        # Step 3: Find the minimum required TP for the assumed path loss
        min_tp_required = assumed_path_loss_db + sensitivity[sf]
        
        # Find the first available TP that meets the requirement
        optimal_tp_for_sf = None
        for tp in tps_dbm_available:
            if tp >= min_tp_required:
                optimal_tp_for_sf = tp
                break
        
        # If a valid TP is found for this SF, calculate energy
        if optimal_tp_for_sf is not None:
            # Calculate Time on Air (ToA)
            de = 1 if sf in [11, 12] else 0  # Low Data Rate Optimization
            t_sym = (2**sf) / bw
            t_preamble = (preamble_symbols + 4.25) * t_sym
            
            payload_symbol_numerator = 8 * payload_bytes - 4 * sf + 28 + (16 if crc_enabled else 0) - (20 if not header_enabled else 0)
            payload_symbol_denominator = 4 * (sf - 2 * de) if (sf - 2*de) != 0 else float('inf')
            cr_multiplier = (coding_rate_val + 4)
            n_payload_sym = 8 + max(math.ceil(payload_symbol_numerator / payload_symbol_denominator) * cr_multiplier, 0)
            t_payload = n_payload_sym * t_sym
            toa_s = t_preamble + t_payload

            # Step 4: Calculate Energy
            tp_watts = 10**((optimal_tp_for_sf - 30) / 10)
            energy_j = tp_watts * toa_s
            
            valid_configs.append({
                "sf": sf,
                "tp": optimal_tp_for_sf,
                "toa_s": toa_s,
                "energy_j": energy_j,
                "tp_watts": tp_watts
            })

    # Step 5: Find the best configuration (minimum energy)
    if not valid_configs:
        print(f"No available (SF, TP) configuration can meet the required path loss of {assumed_path_loss_db} dB.")
        return

    best_config = min(valid_configs, key=lambda x: x['energy_j'])

    # Step 6: Present the final result and the equation
    print("--- Optimal Configuration Analysis ---")
    print(f"For an assumed urban path loss of {assumed_path_loss_db} dB, the most energy-efficient configuration is:")
    print(f"Spreading Factor: {best_config['sf']}")
    print(f"Transmit Power:   {best_config['tp']} dBm\n")
    
    print("This provides the required link reliability with the minimum energy consumption.")
    print("The final energy calculation is:\n")
    
    # Output the final equation with all numbers
    print("Energy (Joules) = Power (Watts) * Time on Air (s)")
    print(f"Energy = 10^(({best_config['tp']} - 30) / 10) * {best_config['toa_s']:.5f}")
    print(f"Energy = {best_config['tp_watts']:.6f} W * {best_config['toa_s']:.5f} s")
    print(f"Energy = {best_config['energy_j']:.6f} Joules ({best_config['energy_j'] * 1000:.2f} mJ)")

# Execute the main function
solve_lorawan_optimization()