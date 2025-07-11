import math

def solve_lorawan_optimization():
    """
    Calculates the optimal LoRaWAN Spreading Factor and Transmit Power to minimize
    energy consumption while meeting a 1% PER requirement in a Rician fading channel.
    """
    # --- Givens and Constants ---
    PAYLOAD_BYTES = 100
    BANDWIDTH_HZ = 125000
    CODING_RATE_CODE = 1  # For CR = 4/5
    HEADER_ENABLED = 1    # LoRaWAN packets have a header
    CRC_ENABLED = 1       # LoRaWAN packets have a CRC
    PREAMBLE_SYMBOLS = 8
    AVAILABLE_TX_POWERS_DBM = [2, 4, 6, 8, 10, 12, 14]

    # --- Derived Constants & Assumptions ---
    # Fading margin for 1% PER (99% PSR) on a Rician channel with K=3dB.
    FADING_MARGIN_DB = 16.51
    # Assumed urban path loss for a meaningful trade-off analysis.
    PATH_LOSS_DB = 115
    # Standard thermal noise at 125kHz BW plus a typical 6dB Gateway Noise Figure.
    RECEIVER_NOISE_FLOOR_DBM = -107

    SF_DEMOD_SNR_DB = {
        7: -7.5,
        8: -10.0,
        9: -12.5,
        10: -15.0,
        11: -17.5,
        12: -20.0
    }

    min_energy = float('inf')
    optimal_config = {}
    
    print("Analyzing all Spreading Factors...")
    print("-" * 50)
    
    # Iterate through all available Spreading Factors
    for sf in range(7, 13):
        # 1. Calculate Time on Air (ToA)
        t_symbol_s = (2**sf) / BANDWIDTH_HZ
        # LowDataRateOptimize flag is enabled if symbol time > 16ms
        de = 1 if t_symbol_s > 0.016 else 0
        
        # Formula for payload symbols
        numerator = 8 * PAYLOAD_BYTES - 4 * sf + 28 + 16 * CRC_ENABLED - 20 * HEADER_ENABLED
        denominator = 4 * (sf - 2 * de)
        
        # Ensure denominator is not zero
        if denominator <= 0:
            continue

        payload_symbols = 8 + max(math.ceil(numerator / denominator) * (CODING_RATE_CODE + 4), 0)
        total_symbols = payload_symbols + PREAMBLE_SYMBOLS + 4.25
        time_on_air_s = total_symbols * t_symbol_s

        # 2. Calculate Required Transmit Power
        snr_demod = SF_DEMOD_SNR_DB[sf]
        required_avg_snr = snr_demod + FADING_MARGIN_DB
        required_tx_power_dbm = required_avg_snr + PATH_LOSS_DB + RECEIVER_NOISE_FLOOR_DBM

        # 3. Find the minimum suitable transmit power from the available list
        chosen_tx_power_dbm = None
        for power in AVAILABLE_TX_POWERS_DBM:
            if power >= required_tx_power_dbm:
                chosen_tx_power_dbm = power
                break
        
        # 4. If a valid power level exists, calculate energy and compare
        if chosen_tx_power_dbm is not None:
            tx_power_mW = 10**(chosen_tx_power_dbm / 10)
            # Energy is proportional to Power * Time on Air
            energy_metric = tx_power_mW * time_on_air_s
            
            if energy_metric < min_energy:
                min_energy = energy_metric
                optimal_config = {
                    "sf": sf,
                    "toa_ms": time_on_air_s * 1000,
                    "snr_demod": snr_demod,
                    "required_avg_snr": required_avg_snr,
                    "required_tx_power_dbm": required_tx_power_dbm,
                    "chosen_tx_power_dbm": chosen_tx_power_dbm,
                    "energy_metric": energy_metric,
                }

    # --- Print Final Result ---
    if not optimal_config:
        print("No valid SF/Power combination found for the assumed path loss.")
    else:
        print("\nOptimal Configuration Found!")
        print("="*50)
        sf = optimal_config['sf']
        chosen_txp = optimal_config['chosen_tx_power_dbm']
        print(f"Optimal Spreading Factor: {sf}")
        print(f"Optimal Transmit Power: {chosen_txp} dBm\n")

        print("Calculation for this optimal setting:")
        # Outputting each number in the final equation as requested
        print(f"1. Base SNR for SF{sf} demodulation: {optimal_config['snr_demod']:.2f} dB")
        print(f"2. Required fading margin for 1% PER (Rician K=3dB): {FADING_MARGIN_DB:.2f} dB")
        print(f"3. Required Average SNR at Receiver = {optimal_config['snr_demod']:.2f} dB + {FADING_MARGIN_DB:.2f} dB = {optimal_config['required_avg_snr']:.2f} dB")
        print(f"4. Required Transmit Power = (Avg SNR) + (Path Loss) + (Noise Floor)")
        print(f"   = {optimal_config['required_avg_snr']:.2f} dB + {PATH_LOSS_DB:.2f} dB + ({RECEIVER_NOISE_FLOOR_DBM:.2f} dBm) = {optimal_config['required_tx_power_dbm']:.2f} dBm")
        print(f"5. Lowest available power that meets requirement: {chosen_txp} dBm")
        chosen_txp_mw = 10**(chosen_txp / 10)
        print(f"6. Time on Air for SF{sf} (100B payload): {optimal_config['toa_ms']:.2f} ms")
        print(f"7. Relative Energy Consumption = {chosen_txp_mw:.2f} mW * {optimal_config['toa_ms']:.2f} ms = {optimal_config['energy_metric']*1000:.2f} (mW*ms)")
        print("="*50)


solve_lorawan_optimization()
<<<10, 10>>>