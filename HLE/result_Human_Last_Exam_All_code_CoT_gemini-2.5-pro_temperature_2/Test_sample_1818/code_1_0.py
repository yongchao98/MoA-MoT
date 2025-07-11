import math

def solve_lorawan_optimization():
    """
    Calculates the optimal LoRaWAN Spreading Factor and Transmit Power
    to minimize energy while ensuring reliability.
    """

    # 1. Given Parameters
    payload_len = 100  # bytes
    bw = 125000  # Hz
    coding_rate_numerator = 4
    coding_rate_denominator = 5
    preamble_symbols = 8
    
    # Available operating points
    spreading_factors = list(range(7, 13))  # SF7 to SF12
    tx_powers_dbm = list(range(2, 15, 2))   # 2dBm to 14dBm in 2dB steps

    # LoRaWAN physical layer parameters for calculations
    # Using accurate ToA values for a 100-byte payload based on online calculators
    # These account for explicit header, CRC, and Low Data Rate Optimization (LDRO)
    toa_s = {
        7: 174.34 / 1000,
        8: 307.71 / 1000,
        9: 553.98 / 1000,
        10: 1026.05 / 1000,
        11: 2052.22 / 1000, # with LDRO
        12: 3776.77 / 1000  # with LDRO
    }

    # 2. Establish Reliability Target (Required SNR)
    # Demodulation sensitivity (bare minimum SNR for reception)
    snr_demod_sensitivity = {
        7: -7.5,
        8: -10.0,
        9: -12.5,
        10: -15.0,
        11: -17.5,
        12: -20.0
    }
    
    # Required link margin for <1% PER in a Rician fading channel
    link_margin_db = 8.0
    
    required_snr = {sf: snr_demod_sensitivity[sf] + link_margin_db for sf in spreading_factors}
    
    # Variables to track the optimal configuration
    optimal_sf = None
    optimal_tx_power = None
    max_efficiency = -1
    
    print("Analyzing all possible (Spreading Factor, Transmit Power) combinations...")
    print("-" * 80)
    print(f"{'SF':>3s} {'TxP (dBm)':>12s} {'ToA (ms)':>12s} {'Energy (mJ)':>15s} {'Link Budget (dB)':>18s} {'Efficiency (LB/E)':>20s}")
    print("-" * 80)

    # 3-5. Iterate through all configurations to find the most efficient one
    for sf in spreading_factors:
        for txp in tx_powers_dbm:
            # Convert dBm to mW for energy calculation
            power_mw = 10**(txp / 10)
            
            # Calculate energy consumption
            energy_mj = power_mw * toa_s[sf] * 1000

            # Calculate link budget
            link_budget = txp - required_snr[sf]
            
            # Calculate efficiency (Link Budget per unit Energy)
            # We only consider configurations with a positive link budget
            if link_budget > 0:
                efficiency = link_budget / energy_mj
            else:
                efficiency = 0
            
            print(f"{sf:3d} {txp:12d} {toa_s[sf]*1000:12.2f} {energy_mj:15.2f} {link_budget:18.2f} {efficiency*1000:20.4f} (x10^-3)")
            
            # Check if this configuration is the new optimum
            if efficiency > max_efficiency:
                max_efficiency = efficiency
                optimal_sf = sf
                optimal_tx_power = txp

    print("-" * 80)
    print("\n--- Conclusion ---")
    print(f"To achieve the reliability goal (PER < 1%) in the specified urban environment, a robust link budget is required.")
    print(f"The analysis identifies the configuration that provides the most link budget for the lowest energy cost (maximum efficiency).")
    print(f"\nThe optimal setting with the highest efficiency ratio (Link Budget / Energy) is:")
    print(f"Spreading Factor: {optimal_sf}")
    print(f"Transmission Power: {optimal_tx_power} dBm")
    
if __name__ == '__main__':
    solve_lorawan_optimization()
<<<SF8, 2dBm>>>