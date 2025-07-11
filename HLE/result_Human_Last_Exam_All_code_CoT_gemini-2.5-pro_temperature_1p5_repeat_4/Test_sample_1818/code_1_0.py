import math

def solve_lorawan_optimization():
    """
    Calculates the optimal LoRaWAN Spreading Factor and Transmission Power
    to minimize energy consumption for a given set of parameters.
    """

    # --- 1. Define Given Parameters ---
    PAYLOAD_BYTES = 100
    BANDWIDTH_HZ = 125000
    CODING_RATE_DENOMINATOR = 5 # From CR = 4/5
    PREAMBLE_LEN = 8
    HEADER_ENABLED = True
    CRC_ON = True
    
    # Available operational parameters
    SPREADING_FACTORS = list(range(7, 13))
    TX_POWERS_DBM = list(range(2, 15, 2))

    # Receiver and Channel characteristics
    # Standard sensitivity (minimum required SNR) for LoRa in AWGN
    SNR_REQ_AWGN_DB = {
        7: -7.5, 8: -10.0, 9: -12.5, 10: -15.0, 11: -17.5, 12: -20.0
    }
    # Fading margin for Rician channel (K=3dB) to achieve PER <= 1%
    FADING_MARGIN_DB = 2.5
    # Assumed gateway receiver noise figure
    NOISE_FIGURE_DB = 6
    # Assumed representative path loss for an urban environment
    PATH_LOSS_DB = 136.0

    # --- 2. Helper Functions ---
    def calculate_toa_s(sf):
        """Calculates the Time on Air (ToA) in seconds for a given SF."""
        # Symbol duration
        t_sym_ms = (2**sf) / (BANDWIDTH_HZ / 1000)
        
        # Low Data Rate Optimization (LDRO) is enabled if symbol time > 16ms
        ldro_enabled = 1 if t_sym_ms > 16 else 0

        # Number of payload symbols
        # Formula from LoRaWAN specification
        numerator = 8 * PAYLOAD_BYTES - 4 * sf + 28 + (16 if CRC_ON else 0) - (20 if HEADER_ENABLED else 0)
        denominator = 4 * (sf - 2 * ldro_enabled)
        
        # Ceil division and multiply by coding rate factor
        payload_symbol_count = 8 + max(0, math.ceil(numerator / denominator) * CODING_RATE_DENOMINATOR)

        # Total symbols (preamble + payload)
        total_symbols = PREAMBLE_LEN + 4.25 + payload_symbol_count
        
        return (total_symbols * t_sym_ms) / 1000 # Convert ms to seconds

    # --- 3. Main Calculation Logic ---
    
    # Calculate receiver noise floor
    noise_floor_dbm = -174 + 10 * math.log10(BANDWIDTH_HZ) + NOISE_FIGURE_DB

    min_energy_joules = float('inf')
    optimal_setting = {}

    print(f"Analyzing for a Path Loss of {PATH_LOSS_DB} dB...")
    print("-" * 50)
    print(f"{'SF':<5} | {'Req. TxP (dBm)':<17} | {'Energy (mJ)':<15}")
    print("-" * 50)
    
    for sf in SPREADING_FACTORS:
        # Calculate the required SNR in the specified Rician channel
        snr_req_rician_db = SNR_REQ_AWGN_DB[sf] + FADING_MARGIN_DB
        
        # Calculate the minimum power required at the receiver
        req_rx_power_dbm = snr_req_rician_db + noise_floor_dbm

        # Calculate the required transmit power to overcome the path loss
        req_tx_power_dbm = req_rx_power_dbm + PATH_LOSS_DB

        # Find the lowest available Tx Power that meets the requirement
        best_tx_power_dbm = None
        for tx_p in TX_POWERS_DBM:
            if tx_p >= req_tx_power_dbm:
                best_tx_power_dbm = tx_p
                break
        
        # If no available power setting works for this SF, skip it
        if best_tx_power_dbm is None:
            print(f"SF{sf:<2} | {'(unachievable)':<17} | {'N/A':<15}")
            continue

        # Calculate energy consumption for this valid (SF, TxP) pair
        tx_power_watts = 10**((best_tx_power_dbm - 30) / 10)
        toa_s = calculate_toa_s(sf)
        energy_joules = tx_power_watts * toa_s

        print(f"SF{sf:<2} | {best_tx_power_dbm:<17.1f} | {energy_joules * 1000:<15.2f}")
        
        # Check if this setting is the new optimum
        if energy_joules < min_energy_joules:
            min_energy_joules = energy_joules
            optimal_setting = {
                "sf": sf,
                "tx_power_dbm": best_tx_power_dbm,
                "tx_power_watts": tx_power_watts,
                "toa_s": toa_s,
                "energy_j": energy_joules
            }

    # --- 4. Output the Final Result ---
    print("-" * 50)
    if not optimal_setting:
        print("\nNo suitable SF and Tx Power combination found for the given path loss.")
    else:
        print("\nOptimal Configuration Found:")
        print(f"  Spreading Factor: SF{optimal_setting['sf']}")
        print(f"  Transmission Power: {optimal_setting['tx_power_dbm']} dBm")
        print(f"  Minimum Energy: {optimal_setting['energy_j'] * 1000:.2f} mJ")
        print("\nFinal Energy Calculation:")
        
        sf_final = optimal_setting['sf']
        txp_final_w = optimal_setting['tx_power_watts']
        toa_final_s = optimal_setting['toa_s']
        eng_final_j = optimal_setting['energy_j']
        
        print(f"  Energy (J) = Transmit Power (W) * Time on Air (s)")
        print(f"  {eng_final_j:.6f} J = {txp_final_w:.6f} W * {toa_final_s:.6f} s")


solve_lorawan_optimization()