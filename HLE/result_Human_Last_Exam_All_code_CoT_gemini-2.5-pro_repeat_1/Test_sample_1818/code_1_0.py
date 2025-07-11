import math

def solve_lorawan_optimization():
    """
    Calculates the optimal Spreading Factor and Transmission Power for a LoRaWAN device
    to minimize energy consumption.
    """
    # Step 1: Define System Parameters
    PL = 100  # Payload in bytes
    BW = 125000  # Bandwidth in Hz
    # For a Coding Rate (CR) of 4/5, the corresponding code (CR_code) is 1.
    CR_code = 1
    H = 1  # Header is enabled
    CRC = 1  # CRC is enabled
    n_preamble = 8  # Number of preamble symbols

    # Define the ranges for SF and Tx Power to be tested
    spreading_factors = list(range(7, 13))  # SF7 to SF12
    tx_powers_dBm = list(range(2, 15, 2))  # 2, 4, ..., 14 dBm

    # Initialize variables to store the optimal results
    min_energy = float('inf')
    optimal_sf = None
    optimal_tx_power_dBm = None
    optimal_toa_s = None
    optimal_power_W = None

    print("Analyzing all combinations of Spreading Factor and Tx Power to find the minimum energy configuration...\n")

    # Iterate through all SF and Tx Power combinations
    for sf in spreading_factors:
        # Step 2: Calculate Time on Air (ToA) for the current SF
        
        # Symbol time in seconds
        t_sym = (2**sf) / BW

        # Determine if Low Data Rate Optimization (DE) is active
        # DE is enabled for SF11 and SF12 at 125kHz BW as T_sym > 16ms
        de = 1 if t_sym > 0.016 else 0

        # Calculate the number of payload symbols using the formula from the LoRa modem specification
        # The term (CR_code + 4) represents the coding rate factor (e.g., 5 for CR 4/5)
        numerator = 8 * PL - 4 * sf + 28 + 16 * CRC - 20 * H
        denominator = 4 * (sf - 2 * de)
        
        # Ensure denominator is not zero
        if denominator == 0:
            continue

        payload_symbols = 8 + math.ceil(max(0, numerator / denominator)) * (CR_code + 4)
        
        # Calculate total Time on Air in seconds
        toa_s = (n_preamble + 4.25 + payload_symbols) * t_sym
        
        for tx_p_dBm in tx_powers_dBm:
            # Step 3: Calculate Energy Consumption
            
            # Convert Tx Power from dBm to Watts
            tx_p_W = (10**(tx_p_dBm / 10)) / 1000
            
            # Calculate energy in Joules
            energy_J = tx_p_W * toa_s
            
            # Step 4: Identify the Optimal Configuration by finding the minimum energy
            if energy_J < min_energy:
                min_energy = energy_J
                optimal_sf = sf
                optimal_tx_power_dBm = tx_p_dBm
                optimal_toa_s = toa_s
                optimal_power_W = tx_p_W

    # Print the final results
    print("--- Optimal Configuration Found ---")
    print(f"Optimal Spreading Factor: SF{optimal_sf}")
    print(f"Optimal Transmission Power: {optimal_tx_power_dBm} dBm")
    print("\nThis combination provides the lowest possible energy consumption. An ADR-enabled network server would assign this setting to a device if the signal strength is sufficient.\n")
    
    print("--- Final Energy Calculation for the Optimal Setting ---")
    print(f"Time on Air for SF{optimal_sf} (ToA): {optimal_toa_s:.6f} seconds")
    print(f"Power for {optimal_tx_power_dBm} dBm (P): {optimal_power_W:.6f} Watts")
    print(f"Resulting Energy (E = P * ToA):")
    print(f"E = {optimal_power_W:.6f} W * {optimal_toa_s:.6f} s")
    print(f"E = {min_energy:.6f} Joules (or {min_energy*1000:.3f} mJ)")

if __name__ == '__main__':
    solve_lorawan_optimization()