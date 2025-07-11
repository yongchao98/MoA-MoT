import math

def calculate_optimal_lorawan_settings():
    """
    Calculates the optimal LoRaWAN Spreading Factor and Transmission Power
    to minimize energy consumption for the given scenario.
    """

    # --- Given Parameters ---
    PL = 100          # Payload in bytes
    BW = 125000       # Bandwidth in Hz
    CR_val = 4/5      # Coding Rate
    N_preamble = 8    # Preamble symbols
    H = 0             # Header is explicit (H=0)
    CRC = 1           # CRC is on (1)

    # Coding Rate code for the formula (1 for 4/5, 2 for 4/6, etc.)
    # CR = 4 / (4 + cr_code)
    cr_code = int(4 / CR_val - 4)

    # SNR demodulation thresholds in dB for each Spreading Factor
    SNR_THRESHOLDS = {
        7: -7.5,
        8: -10.0,
        9: -12.5,
        10: -15.0,
        11: -17.5,
        12: -20.0
    }
    
    # Available transmit powers in dBm
    TX_POWERS_DBM = list(range(2, 15, 2))

    print("Plan: Calculate the 'Energy Factor' for each Spreading Factor (SF).")
    print("The factor is proportional to the energy required to transmit one packet.")
    print("Energy Factor = 10^(SNR_threshold / 10) * Time_on_Air\n")
    print("Lower is more energy-efficient.\n")
    
    results = {}

    # --- Iterate through each Spreading Factor ---
    for sf in range(7, 13):
        # --- 1. Calculate Time on Air (ToA) ---
        
        # Low Data Rate Optimization (DE) is enabled for SF11 and SF12 on 125kHz
        de = 1 if sf >= 11 else 0
        
        # Symbol duration
        t_sym = (2**sf) / BW
        
        # Calculate number of payload symbols
        payload_sym_num_numerator = 8 * PL - 4 * sf + 28 + 16 * CRC - 20 * H
        payload_sym_num_denominator = 4 * (sf - 2 * de)
        
        # Note: The expression must be an integer, hence math.ceil
        payload_sym_num = 8 + max(0, math.ceil(payload_sym_num_numerator / payload_sym_num_denominator)) * (cr_code + 4)

        # Total time on air in seconds
        toa = (N_preamble + 4.25 + payload_sym_num) * t_sym
        
        # --- 2. Calculate Energy Factor ---
        snr_threshold_db = SNR_THRESHOLDS[sf]
        # Convert SNR from dB to linear scale
        gamma_threshold_linear = 10**(snr_threshold_db / 10)
        
        # The energy factor is the product of the linear SNR and ToA
        energy_factor = gamma_threshold_linear * toa
        
        results[sf] = {'toa': toa, 'factor': energy_factor}
        
        print(f"--- SF{sf} ---")
        print(f"Time on Air (ToA)       = {toa:.4f} seconds")
        print(f"Required SNR Threshold  = {snr_threshold_db} dB")
        # Explicitly show the numbers used in the final equation
        print(f"Energy Factor           = 10^({snr_threshold_db} / 10) * {toa:.4f} s = {energy_factor:.5f}\n")
    
    # --- 3. Determine the optimal SF ---
    # Find the SF with the minimum energy factor
    optimal_sf = min(results, key=lambda sf: results[sf]['factor'])
    
    # --- 4. Determine the optimal TxP ---
    # To minimize energy, the ADR should select the lowest possible power
    # that satisfies the link budget. We select the minimum available power.
    optimal_txp = min(TX_POWERS_DBM)

    print("--- Conclusion ---")
    print(f"The most energy-efficient Spreading Factor is SF{optimal_sf}, as it has the lowest Energy Factor.")
    print(f"To minimize energy consumption, the network server should use this SF and select the lowest possible")
    print(f"transmit power that ensures a PER <= 1%. The minimum available power is {optimal_txp} dBm.")
    print("\nOptimal Configuration:")
    print(f"Spreading Factor: {optimal_sf}")
    print(f"Transmission Power: {optimal_txp} dBm")


if __name__ == '__main__':
    calculate_optimal_lorawan_settings()
<<<Spreading Factor: 8, Transmission Power: 2 dBm>>>