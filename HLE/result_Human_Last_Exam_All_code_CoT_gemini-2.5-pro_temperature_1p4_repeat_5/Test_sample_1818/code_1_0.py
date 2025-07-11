import math

def calculate_energy_optimization():
    """
    Calculates the optimal Spreading Factor and Transmit Power for a LoRaWAN device
    to minimize energy while meeting a reliability target in an urban environment.
    """
    # 1. Define constants from the problem
    transmit_powers_dbm = [2, 4, 6, 8, 10, 12, 14]
    spreading_factors = [7, 8, 9, 10, 11, 12]
    bw = 125000  # Bandwidth in Hz
    payload_bytes = 100
    coding_rate_val = 4/5
    codr = 5 # Denominator of the coding rate 4/5
    preamble_symbols = 8
    header_implicit = 0 # 0 for explicit header
    crc_on = 16 # CRC is on (16 bit)

    # Required SNR for demodulation at low PER for each SF (dB)
    snr_required = {7: -7.5, 8: -10, 9: -12.5, 10: -15, 11: -17.5, 12: -20}

    # 2. Model the link requirements
    # Fading margin for Rician K=3dB at 99% reliability (for 1% PER)
    fading_margin_db = 8.2
    
    # Receiver gateway parameters to calculate sensitivity
    thermal_noise_dbm_hz = -174 # dBm/Hz at room temperature
    noise_bw_db = 10 * math.log10(bw)
    gateway_noise_figure_db = 6 # Typical for a LoRa gateway

    # Assumed path loss for an "urban environment"
    assumed_path_loss_db = 135
    
    print(f"Analysis for an assumed path loss of {assumed_path_loss_db} dB.\n")

    candidate_configs = []

    # 3. Find candidate configurations for the assumed path loss
    for sf in spreading_factors:
        # Calculate effective required SNR including fading margin
        effective_snr_req_db = snr_required[sf] + fading_margin_db
        
        # Calculate receiver sensitivity for this SF
        # Sensitivity = ThermalNoise + NoiseBandwidth + NoiseFigure + RequiredSNR
        sensitivity_dbm = thermal_noise_dbm_hz + noise_bw_db + gateway_noise_figure_db + effective_snr_req_db

        # Find the minimum TP to overcome the path loss
        min_tp_needed = None
        for tp in transmit_powers_dbm:
            link_budget = tp - sensitivity_dbm
            if link_budget >= assumed_path_loss_db:
                min_tp_needed = tp
                break # Found the lowest possible power for this SF

        if min_tp_needed is not None:
            candidate_configs.append({'sf': sf, 'tp': min_tp_needed})
        else:
            # This SF cannot meet the path loss requirement even at max power
            pass

    print("Finding optimal configuration from the following candidates:")

    # 4. Calculate Time on Air and Energy for each candidate
    optimal_config = None
    min_energy = float('inf')

    for config in candidate_configs:
        sf = config['sf']
        tp_dbm = config['tp']
        
        # Time on Air (ToA) Calculation based on Semtech formula
        low_data_rate_optimize = 1 if sf in [11, 12] else 0
        
        numerator = 8 * payload_bytes - 4 * sf + 28 + crc_on - 20 * header_implicit
        denominator = 4 * (sf - 2 * low_data_rate_optimize)
        
        # Number of payload symbols
        payload_symbol_count = 8 + max(0, math.ceil(numerator / denominator) * codr)
        
        symbol_time_s = (2**sf) / bw
        preamble_time_s = (preamble_symbols + 4.25) * symbol_time_s
        payload_time_s = payload_symbol_count * symbol_time_s
        toa_s = preamble_time_s + payload_time_s

        # Energy Calculation (relative value)
        # Energy = Power (mW) * Time (s)
        tp_mw = 10**(tp_dbm / 10)
        energy_cost = tp_mw * toa_s

        # Store results for the current candidate
        config['toa_s'] = toa_s
        config['energy_cost'] = energy_cost
        
        print(f"\n- Candidate: SF{sf}, {tp_dbm} dBm")
        print(f"  - Time on Air: {toa_s:.3f} s")
        print(f"  - Energy Equation: {tp_mw:.2f} mW * {toa_s:.3f} s")
        print(f"  - Calculated Energy Cost: {energy_cost:.2f} mJ (relative units)")

        if energy_cost < min_energy:
            min_energy = energy_cost
            optimal_config = config
            
    # 5. Announce the final result
    if optimal_config:
        print("\n------------------------------------------------------------")
        print("Comparison complete.")
        print(f"The optimal configuration with the lowest energy consumption is:")
        print(f"Spreading Factor {optimal_config['sf']} and Transmission Power {optimal_config['tp']} dBm")
        print("This configuration meets the 1% PER requirement for the assumed urban path loss with the minimum energy expenditure.")
        print("------------------------------------------------------------")
        # Returning the final answer in the required format
        # print(f"\n<<<SF{optimal_config['sf']}, {optimal_config['tp']}dBm>>>")
    else:
        print("\nNo suitable configuration found for the assumed path loss.")


calculate_energy_optimization()
# Final answer based on the script's logic and calculation.
# For a path loss of 135dB, the candidates are:
# (SF9, 14dBm): E = 10**(1.4) * 0.554 = 25.12 * 0.554 = 13.92
# (SF10, 12dBm): E = 10**(1.2) * 1.026 = 15.85 * 1.026 = 16.26
# (SF11, 10dBm): E = 10**(1.0) * 2.216 = 10 * 2.216 = 22.16
# (SF12, 8dBm): E = 10**(0.8) * 3.940 = 6.31 * 3.940 = 24.86
# The lowest energy is 13.92, which corresponds to SF9 and 14dBm.
final_sf = 9
final_tp = 14
print(f"Final Answer: The optimal configuration is Spreading Factor {final_sf} and Transmission Power {final_tp} dBm.")