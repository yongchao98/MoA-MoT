import math

def solve_lora_optimization():
    """
    Calculates the optimal LoRaWAN Spreading Factor and Transmit Power to minimize
    energy consumption for a given set of parameters.
    """

    # --- Problem Parameters ---
    PAYLOAD_BYTES = 100
    BANDWIDTH_HZ = 125000
    CODING_RATE = "4/5" # This implies a CR value of 1 in the formula
    CODING_RATE_VAL = 1
    PREAMBLE_LEN = 8
    # For standard LoRaWAN packets, the header is explicit, so Header_Enabled (H) is 0
    HEADER_ENABLED = 0
    # Available transmit power levels in dBm
    TX_POWERS_DBM = list(range(2, 15, 2))

    # Spreading Factors and their typical required SNR (dB) for ~1% PER
    # These are the demodulation thresholds at the gateway.
    SF_SNR_MAP = {
        7: -7.5,
        8: -10.0,
        9: -12.5,
        10: -15.0,
        11: -17.5,
        12: -20.0
    }

    print("Analyzing LoRaWAN Energy Consumption for a 100-byte Payload")
    print("Parameters: Bandwidth={} kHz, Coding Rate={}, Header=Explicit".format(int(BANDWIDTH_HZ / 1000), CODING_RATE))
    print("-" * 75)
    print("{:<5} | {:<10} | {:<12} | {:<15} | {:<20}".format("SF", "SNR (dB)", "ToA (ms)", "Rel. Power Req.", "Relative Energy"))
    print("-" * 75)

    results = []

    # Iterate through each Spreading Factor to find the most energy-efficient one
    for sf, snr_db in SF_SNR_MAP.items():
        # --- 1. Calculate Time on Air (ToA) ---
        t_symbol = (2**sf) / BANDWIDTH_HZ
        
        # Data Rate Optimization is used for SF11 and SF12
        de = 1 if sf >= 11 else 0

        # Calculate the number of payload symbols using the standard formula
        # Numerator: 8*PL - 4*SF + 28 + 16(CRC) - 20*H
        numerator = 8 * PAYLOAD_BYTES - 4 * sf + 28 + 16 - (20 * HEADER_ENABLED)
        denominator = 4 * (sf - 2 * de)
        
        # Symbol count must be an integer
        payload_symbol_count = 8 + math.ceil(max(0, numerator) / denominator) * (CODING_RATE_VAL + 4)
        
        # Total symbols including preamble and header
        total_symbols = PREAMBLE_LEN + 4.25 + payload_symbol_count
        
        # Final ToA in seconds
        toa_s = total_symbols * t_symbol

        # --- 2. Calculate Relative Energy Metric ---
        # The required transmit power is proportional to 10^(SNR/10)
        relative_power_linear = 10**(snr_db / 10)
        
        # The relative energy is proportional to Power * Time
        relative_energy = relative_power_linear * toa_s

        results.append({
            "sf": sf,
            "snr_db": snr_db,
            "toa_ms": toa_s * 1000,
            "relative_power": relative_power_linear,
            "relative_energy": relative_energy,
            "toa_s": toa_s
        })

    # Print the calculated table for analysis
    for res in results:
        print("{:<5} | {:<10.1f} | {:<12.2f} | {:<15.4f} | {:<20.6f}".format(
            res['sf'], res['snr_db'], res['toa_ms'], res['relative_power'], res['relative_energy']
        ))

    # --- 3. Determine the Optimal Configuration ---
    # Find the configuration with the minimum relative energy
    optimal_config = min(results, key=lambda x: x['relative_energy'])
    optimal_sf = optimal_config['sf']

    # For energy minimization, the lowest available transmit power is the target
    optimal_tx_power = min(TX_POWERS_DBM)

    print("-" * 75)
    print("\nConclusion:")
    print(f"The analysis shows that SF{optimal_sf} has the lowest relative energy consumption.")
    print("This means for any given distance, it is the most energy-efficient spreading factor for this payload size.")
    print("\nTo achieve the lowest energy consumption, the network server should assign:")
    print(f"  - Optimal Spreading Factor: SF{optimal_sf}")
    print(f"  - Optimal Transmission Power: {optimal_tx_power} dBm (as the lowest available and sufficient power)")

    # --- 4. Output the final equation with numbers ---
    print("\nFinal Calculation for Optimal Configuration:")
    print("The optimal setting is found by minimizing the relative energy equation:")
    print("E_relative = 10^(SNR_req / 10) * Time_on_Air")
    print("\nFor the optimal SF{}:".format(optimal_sf))
    print("E_relative = 10^({:.1f} dB / 10) * {:.4f} s".format(optimal_config['snr_db'], optimal_config['toa_s']))
    print("E_relative = {:.4f} * {:.4f} s = {:.6f} (arbitrary energy units)".format(
        optimal_config['relative_power'],
        optimal_config['toa_s'],
        optimal_config['relative_energy']
    ))


if __name__ == '__main__':
    solve_lora_optimization()
    print("\n<<<SF8, 2 dBm>>>")