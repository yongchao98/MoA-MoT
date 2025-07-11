import math
import numpy as np
from scipy.special import marcumq
from scipy.optimize import brentq

def calculate_toa_ms(payload_len, sf, bw_khz=125, cr_code=1, crc=True, header=True, preamble_len=8):
    """
    Calculates the Time on Air (ToA) for a LoRa packet in milliseconds.
    """
    bw_hz = bw_khz * 1000
    # Low Data Rate Optimization for SF11 and SF12
    de = 1 if sf >= 11 else 0
    # Symbol time in ms
    t_sym = (2**sf / bw_hz) * 1000

    # Header is enabled (H=0), CRC is on (16)
    h = 0 if header else 1
    crc_val = 16 if crc else 0

    # Numerator of the payload symbols calculation
    payload_num = 8 * payload_len - 4 * sf + 28 + crc_val - 20 * h
    # Denominator of the payload symbols calculation
    payload_den = 4 * (sf - 2 * de)

    # Calculate number of payload symbols
    if payload_num <= 0:
        n_payload = 0
    else:
        n_payload = (cr_code + 4) * math.ceil(payload_num / payload_den)

    n_total_symbols = preamble_len + 4.25 + 8 + n_payload
    toa = n_total_symbols * t_sym
    return toa

def find_assigned_tx_power(required_power, available_powers):
    """
    Finds the lowest available transmit power that meets the requirement.
    """
    for power in available_powers:
        if power >= required_power:
            return power
    # If no power level is sufficient
    return float('inf')

# --- 1. Define Parameters ---
PAYLOAD_BYTES = 100
BANDWIDTH_KHZ = 125
CODING_RATE_CODE = 1  # For CR = 4/5
K_FACTOR_DB = 3
TARGET_PER = 0.01
TARGET_PSR = 1.0 - TARGET_PER

# Environment and System Parameters
ASSUMED_PATH_LOSS_DB = 120 # Typical value for urban NLOS
GATEWAY_NOISE_FIGURE_DB = 6
NOISE_FLOOR_DBM = -174 + 10 * np.log10(BANDWIDTH_KHZ * 1000) + GATEWAY_NOISE_FIGURE_DB

# Available device settings
SPREADING_FACTORS = list(range(7, 13))
AVAILABLE_TX_POWERS_DBM = [2, 4, 6, 8, 10, 12, 14]

# SNR thresholds for each SF (typical values for 125kHz BW)
SNR_THRESHOLDS_DB = {
    7: -7.5,
    8: -10.0,
    9: -12.5,
    10: -15.0,
    11: -17.5,
    12: -20.0,
}

# --- 2. Calculate Required Fading Margin for Rician Channel ---
K_lin = 10**(K_FACTOR_DB / 10)
# Parameter 'a' for Marcum Q-function
a_marcum = np.sqrt(2 * K_lin)
# Function to find the root for 'b'
func_to_solve = lambda b: marcumq(a_marcum, b, 1) - TARGET_PSR
# Solve for 'b' where the function equals zero
b_solved = brentq(func_to_solve, 0, a_marcum * 2)
# Calculate the required linear margin over the threshold
margin_factor = (2 * (K_lin + 1)) / (b_solved**2)
FADING_MARGIN_DB = 10 * np.log10(margin_factor)

print("--- System Analysis ---")
print(f"Assumed Path Loss: {ASSUMED_PATH_LOSS_DB} dB")
print(f"Gateway Noise Floor: {NOISE_FLOOR_DBM:.1f} dBm")
print(f"Rician Fading Margin for 99% PSR (K={K_FACTOR_DB}dB): {FADING_MARGIN_DB:.1f} dB\n")

# --- 3. Iterate and Find Optimal Configuration ---
results = []
print("--- Energy Calculation for Each Spreading Factor ---")
for sf in SPREADING_FACTORS:
    # Calculate Time on Air
    toa_ms = calculate_toa_ms(PAYLOAD_BYTES, sf)

    # Calculate required average SNR at the receiver
    snr_thresh_db = SNR_THRESHOLDS_DB[sf]
    required_snr_avg_db = snr_thresh_db + FADING_MARGIN_DB

    # Calculate required transmit power
    required_tx_power_dbm = ASSUMED_PATH_LOSS_DB + NOISE_FLOOR_DBM + required_snr_avg_db

    # Find the discrete power level to be assigned by the server
    assigned_tx_power_dbm = find_assigned_tx_power(required_tx_power_dbm, AVAILABLE_TX_POWERS_DBM)

    # Calculate energy metric (ToA * Power_linear)
    if assigned_tx_power_dbm == float('inf'):
        energy_metric = float('inf')
        print(f"SF{sf}: Required TxP ({required_tx_power_dbm:.1f} dBm) exceeds maximum available power (14 dBm).")
    else:
        assigned_tx_power_mw = 10**(assigned_tx_power_dbm / 10)
        energy_metric = toa_ms * assigned_tx_power_mw
        # Output each number in the final equation
        print(f"SF{sf}: Required TxP={required_tx_power_dbm:.1f} dBm -> Assigned TxP={assigned_tx_power_dbm} dBm. "
              f"Energy = {toa_ms:.1f} ms * {assigned_tx_power_mw:.2f} mW = {energy_metric:.0f}")

    results.append({
        'sf': sf,
        'toa_ms': toa_ms,
        'required_tx_power_dbm': required_tx_power_dbm,
        'assigned_tx_power_dbm': assigned_tx_power_dbm,
        'energy_metric': energy_metric
    })

# --- 4. Determine and Print the Optimal Result ---
optimal_result = min(results, key=lambda x: x['energy_metric'])

print("\n--- Optimal Configuration ---")
if optimal_result['energy_metric'] == float('inf'):
    print("No configuration can meet the PER requirement under the assumed path loss.")
else:
    optimal_sf = optimal_result['sf']
    optimal_tx_p = optimal_result['assigned_tx_power_dbm']
    print(f"The most energy-efficient configuration to meet the 1% PER requirement is:")
    print(f"Spreading Factor: SF{optimal_sf}")
    print(f"Transmission Power: {optimal_tx_p} dBm")
    
    # Final answer in the specified format
    # print(f"\n<<<{optimal_sf}, {optimal_tx_p} dBm>>>")
