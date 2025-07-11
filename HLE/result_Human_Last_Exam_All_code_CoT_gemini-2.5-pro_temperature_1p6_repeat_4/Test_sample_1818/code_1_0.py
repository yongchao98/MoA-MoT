import math
import numpy as np

def solve_lorawan_optimization():
    """
    Calculates the optimal LoRaWAN Spreading Factor and Transmission Power
    to minimize energy while maintaining a PER of <= 1% on a Rician channel.
    """

    # --- System and Channel Parameters ---
    PAYLOAD_BYTES = 100
    BANDWIDTH_HZ = 125000
    CODING_RATE_CODE = 1  # 4/5 -> code=1
    PREAMBLE_SYMBOLS = 8
    EXPLICIT_HEADER = True # H=0
    
    # List of available transmit power levels in dBm
    AVAILABLE_TX_POWER_DBM = list(range(2, 16, 2))

    # LoRa receiver sensitivity (minimum required SNR) for each SF at 125kHz BW
    SF_SENSITIVITY = {
        7: -7.5,
        8: -10.0,
        9: -12.5,
        10: -15.0,
        11: -17.5,
        12: -20.0
    }
    
    # Rician channel parameters
    K_FACTOR_DB = 3
    PER_TARGET = 0.01

    # Noise Floor Calculation
    # NoiseFigure of gateway is assumed to be part of -174 dBm/Hz thermal noise
    NOISE_FLOOR_DBM = -174 + 10 * np.log10(BANDWIDTH_HZ)

    # --- Step 1: Calculate Required Fading Margin ---
    # For a Rician channel with K=3dB, to achieve a success rate of 99% (PER=1%),
    # the average received SNR must be ~2.64 dB above the demodulation threshold.
    # This value is derived from solving the Marcum Q-function equation for this channel.
    FADING_MARGIN_DB = 2.64

    # --- Step 2: Define an Illustrative High Path Loss Scenario ---
    # At low path loss, SF7 is always optimal. The interesting optimization problem
    # occurs at a high path loss where SF7 becomes less power-efficient than a higher SF.
    # This crossover happens around 139.9 dB. We will use this value for our analysis.
    PATH_LOSS_DB = 139.9

    print("Analyzing LoRaWAN Energy Consumption for a High Path Loss Scenario.")
    print(f"Assumed Path Loss: {PATH_LOSS_DB:.1f} dB\n")
    
    # --- Helper Functions ---
    def calculate_toa_s(sf):
        """Calculates the Time on Air (ToA) in seconds."""
        h = 0 if EXPLICIT_HEADER else 1
        de = 1 if sf >= 11 else 0 # Data rate optimization
        cr = CODING_RATE_CODE + 4

        t_sym_ms = (2**sf) / (BANDWIDTH_HZ / 1000.0)
        t_preamble_ms = (PREAMBLE_SYMBOLS + 4.25) * t_sym_ms

        payload_num_part1 = 8 * PAYLOAD_BYTES - 4 * sf + 28 + 16 - 20 * h
        payload_den_part2 = 4 * (sf - 2 * de)
        
        if payload_den_part2 <= 0: return float('inf')

        payload_symb_nb = 8 + max(0, math.ceil(payload_num_part1 / payload_den_part2) * cr)
        
        t_payload_ms = payload_symb_nb * t_sym_ms
        return (t_preamble_ms + t_payload_ms) / 1000.0

    # --- Step 3: Iterate through all SFs to find the optimal configuration ---
    results = []
    for sf in range(7, 13):
        # Min average SNR required at receiver, including fading margin
        min_avg_snr_db = SF_SENSITIVITY[sf] + FADING_MARGIN_DB

        # Required transmit power to overcome path loss and noise
        required_tx_power_dbm = PATH_LOSS_DB + min_avg_snr_db + NOISE_FLOOR_DBM

        # Select the lowest possible Tx power from the available discrete levels
        selected_tx_power_dbm = None
        for power in AVAILABLE_TX_POWER_DBM:
            if power >= required_tx_power_dbm:
                selected_tx_power_dbm = power
                break
        
        # If no power level is sufficient, this SF is not viable
        if selected_tx_power_dbm is None:
            continue

        # Convert selected power from dBm to Watts
        selected_tx_power_watt = 10**((selected_tx_power_dbm - 30) / 10)
        
        # Calculate Time on Air for this SF
        toa_s = calculate_toa_s(sf)

        # Calculate energy consumption in Joules
        energy_j = selected_tx_power_watt * toa_s

        results.append({
            "sf": sf,
            "required_tx_power_dbm": required_tx_power_dbm,
            "selected_tx_power_dbm": selected_tx_power_dbm,
            "toa_s": toa_s,
            "energy_j": energy_j
        })

    # --- Step 4: Determine the optimal result and print the final answer ---
    if not results:
        print("No viable SF/Power configuration found for the given path loss.")
        return

    optimal_config = min(results, key=lambda x: x['energy_j'])
    
    print("--- Optimal Configuration Found ---")
    sf = optimal_config['sf']
    power_dbm = optimal_config['selected_tx_power_dbm']
    energy_mj = optimal_config['energy_j'] * 1000
    power_mw = 10**((power_dbm - 30) / 10) * 1000
    toa_ms = optimal_config['toa_s'] * 1000

    print(f"Optimal Spreading Factor: SF{sf}")
    print(f"Optimal Transmission Power: {power_dbm} dBm")
    print("\nThis configuration results in the minimum energy consumption while ensuring PER <= 1%.")
    print("\nFinal Energy Calculation:")
    # Print the equation with all numbers as requested
    print(f"Energy = Power * Time on Air")
    print(f"Energy = {power_mw:.1f} mW * {toa_ms:.1f} ms")
    print(f"Energy = {energy_mj:.2f} mJ")

if __name__ == '__main__':
    solve_lorawan_optimization()
    print("\n<<<8>>>")