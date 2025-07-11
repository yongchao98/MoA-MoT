import math

def calculate_energy_and_find_optimum():
    """
    Calculates the optimal SF and TxP for a LoRaWAN device based on given parameters.
    """
    # --- Given Parameters ---
    payload_len = 100  # bytes
    tx_powers_dbm = [2, 4, 6, 8, 10, 12, 14]
    spreading_factors = [7, 8, 9, 10, 11, 12]
    bandwidth = 125000  # Hz
    coding_rate_val = 4/5
    preamble_len = 8
    header_enabled = True # H=0 in formula

    # --- Assumptions ---
    # Fading margin for Rician K=3dB at PER 1% (engineering estimate)
    fading_margin_db = 4
    # Assumed path loss for an "urban environment"
    assumed_path_loss_db = 143
    # Receiver noise floor for 125kHz BW and 3dB Noise Figure
    noise_floor_dbm = -120

    # LoRaWAN receiver sensitivity (SNR) for 125kHz BW (standard values)
    snr_sensitivity = {
        7: -7.5, 8: -10, 9: -12.5, 10: -15, 11: -17.5, 12: -20
    }

    # --- Calculation ---
    
    valid_configs = []

    print("Analyzing all SF and Tx Power combinations...\n")

    for sf in spreading_factors:
        # Calculate required SNR at the receiver
        required_snr = snr_sensitivity[sf] + fading_margin_db
        
        # Calculate Time on Air (ToA)
        # Using the formula from the LoRaWAN specification
        h = 0 if header_enabled else 1
        de = 1 if sf >= 11 else 0 # Data Rate Optimization
        cr_code = 1 # For CR = 4/5
        
        payload_symbol_count_numerator = 8 * payload_len - 4 * sf + 28 + 16 - 20 * h
        payload_symbol_count_denominator = 4 * (sf - 2 * de)
        
        payload_symbol_count = 8 + max(math.ceil(payload_symbol_count_numerator / payload_symbol_count_denominator) * (cr_code + 4), 0)

        symbol_duration = (2**sf) / bandwidth
        preamble_duration = (preamble_len + 4.25) * symbol_duration
        payload_duration = payload_symbol_count * symbol_duration
        time_on_air = preamble_duration + payload_duration

        for txp in tx_powers_dbm:
            # Maximum path loss this configuration can tolerate
            max_pl_supported = txp - (required_snr + noise_floor_dbm)

            # Check if this configuration is viable for the assumed environment
            if max_pl_supported >= assumed_path_loss_db:
                # Calculate energy consumption in Joules
                txp_watt = 10**((txp - 30) / 10)
                energy_joules = time_on_air * txp_watt
                valid_configs.append({
                    "sf": sf,
                    "txp": txp,
                    "energy": energy_joules,
                    "toa": time_on_air,
                    "power_w": txp_watt
                })

    if not valid_configs:
        print(f"No valid configuration found for the assumed path loss of {assumed_path_loss_db} dB.")
        return

    # Find the configuration with the minimum energy
    optimal_config = min(valid_configs, key=lambda x: x['energy'])

    # --- Output the result ---
    print("--- Optimal Configuration Found ---")
    print(f"Optimal Spreading Factor (SF): {optimal_config['sf']}")
    print(f"Optimal Transmit Power (TxP): {optimal_config['txp']} dBm")
    print("\nThis configuration is the most energy-efficient choice that overcomes the")
    print(f"assumed {assumed_path_loss_db} dB path loss while meeting the 1% PER requirement.\n")

    print("Energy Calculation for the Optimal Configuration:")
    # Using 'g' for general-purpose number formatting to avoid unnecessary trailing zeros
    print(f"Energy (J) = Time on Air (s) * Transmit Power (W)")
    print(f"Lowest Energy = {optimal_config['toa']:.3f} s * {optimal_config['power_w']:.4f} W = {optimal_config['energy']:.4f} J")


if __name__ == '__main__':
    calculate_energy_and_find_optimum()
<<<SF10, 12 dBm>>>