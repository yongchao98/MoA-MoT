import math

def solve_lorawan_optimization():
    """
    Calculates the optimal LoRaWAN Spreading Factor and Transmit Power
    to minimize energy consumption for a given scenario.
    """

    # --- Given Parameters ---
    PAYLOAD_BYTES = 100
    BANDWIDTH_HZ = 125000
    CODING_RATE_STR = "4/5" # Corresponds to CR code 1 in the formula
    PREAMBLE_SYMBOLS = 8
    HAS_HEADER = True
    HAS_CRC = True
    TX_POWERS_DBM = list(range(2, 15, 2)) # 2, 4, 6, 8, 10, 12, 14 dBm

    # --- Channel and System Parameters ---
    # For Rician channel with K=3dB and PER<=1%, a fading margin is needed.
    FADING_MARGIN_DB = 4.5
    # Thermal noise floor for 125kHz BW = -174 + 10*log10(BW)
    NOISE_FLOOR_DBM = -174 + 10 * math.log10(BANDWIDTH_HZ)
    # Assumed typical urban path loss
    PATH_LOSS_DB = 137.0

    # SF-specific data: (SF, SNR_threshold_AWGN, LowDataRateOptimize_flag)
    SF_DATA = {
        7: (-7.5, False),
        8: (-10.0, False),
        9: (-12.5, False),
        10: (-15.0, False),
        11: (-17.5, True),
        12: (-20.0, True),
    }

    print("--- LoRaWAN Energy Optimization Analysis ---")
    print(f"Assuming a fixed path loss of {PATH_LOSS_DB} dB.\n")
    
    results = []

    for sf, (snr_awgn, de_flag) in SF_DATA.items():
        # 1. Calculate Time on Air (ToA)
        tsym_s = (2**sf) / BANDWIDTH_HZ
        de = 1 if de_flag else 0
        
        # Formula for number of payload symbols
        payload_numerator = 8 * PAYLOAD_BYTES - 4 * sf + 28 + (16 if HAS_CRC else 0) - (20 if HAS_HEADER else 0)
        payload_denominator = 4 * (sf - 2 * de)
        
        # For CR 4/5, the multiplier is 5 (CR_code 1 + 4)
        n_payload_sym = 8 + max(0, math.ceil(payload_numerator / payload_denominator) * 5)
        
        # Total symbols = Preamble (8) + SFD (4.25) + Payload
        n_total_sym = n_payload_sym + PREAMBLE_SYMBOLS + 4.25
        toa_s = n_total_sym * tsym_s

        # 2. Calculate Required Receiver Sensitivity
        snr_required_db = snr_awgn - FADING_MARGIN_DB # Semtech SNR is negative
        sensitivity_dbm = NOISE_FLOOR_DBM + snr_required_db

        # 3. Determine Required Transmit Power for the assumed path loss
        required_tx_power_dbm = PATH_LOSS_DB + sensitivity_dbm
        
        # Find the smallest available power level that meets the requirement
        optimal_tx_power_dbm = None
        for pwr in TX_POWERS_DBM:
            if pwr >= required_tx_power_dbm:
                optimal_tx_power_dbm = pwr
                break
        
        # 4. Calculate Energy if a viable power level was found
        if optimal_tx_power_dbm is not None:
            tx_power_mw = 10**(optimal_tx_power_dbm / 10)
            energy_mj = tx_power_mw * toa_s * 1000 # Convert Joules to milliJoules
            
            results.append({
                "sf": sf,
                "toa_ms": toa_s * 1000,
                "sensitivity_dbm": sensitivity_dbm,
                "required_tx_dbm": required_tx_power_dbm,
                "chosen_tx_dbm": optimal_tx_power_dbm,
                "energy_mj": energy_mj
            })

    # 5. Analyze results and find the minimum energy configuration
    print("--- Viable Configurations Analysis ---")
    if not results:
        print("No viable configuration found for the given path loss.")
        return

    print(f"{'SF':<4} | {'ToA (ms)':>10} | {'Sensitivity (dBm)':>18} | {'Required Power (dBm)':>22} | {'Chosen Power (dBm)':>20} | {'Energy (mJ)':>13}")
    print("-" * 115)
    for res in results:
        print(f"{res['sf']:<4} | {res['toa_ms']:10.2f} | {res['sensitivity_dbm']:18.2f} | {res['required_tx_dbm']:22.2f} | {res['chosen_tx_dbm']:>20} | {res['energy_mj']:13.2f}")

    min_energy_config = min(results, key=lambda x: x['energy_mj'])

    print("\n--- Optimal Configuration ---")
    print("The configuration that meets the PER requirement with the minimum energy consumption is:")
    
    # Final Answer Output
    optimal_sf = min_energy_config['sf']
    optimal_power = min_energy_config['chosen_tx_dbm']
    print(f"\nSpreading Factor (SF): {optimal_sf}")
    print(f"Transmission Power: {optimal_power} dBm")
    
    final_equation = f"Energy({min_energy_config['energy_mj']:.2f} mJ) = Power({10**(optimal_power/10):.2f} mW) * Time({min_energy_config['toa_ms']/1000:.3f} s)"
    print(f"\nFinal Equation: {final_equation}")
    print(f"To satisfy the output format, the numbers in the final equation are: {min_energy_config['energy_mj']:.2f}, {10**(optimal_power/10):.2f}, and {min_energy_config['toa_ms']/1000:.3f}")
    
# Execute the solver
solve_lorawan_optimization()
