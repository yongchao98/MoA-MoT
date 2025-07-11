import math

def db_to_linear(db_value):
    """Converts a dB value to a linear ratio."""
    return 10**(db_value / 10)

def linear_to_db(linear_value):
    """Converts a linear ratio to a dB value."""
    return 10 * math.log10(linear_value)

# --- Given Parameters ---
# Transmitter (Tx)
tx_power_dbm = 30.0
tx_antenna_gain_db = 20.0
tx_antenna_loss_db = 1.0
tx_filter_loss_db = 1.0
tx_cable_loss_db = 1.0

# Path
frequency_ghz = 24.0
distance_km = 10.0

# Receiver (Rx)
rx_antenna_gain_db = 1.0
rx_antenna_loss_db = 0.5
rx_input_filter_loss_db = 1.0
lna_gain_db = 36.0
lna_nf_db = 2.0
mixer_loss_db = 9.0
if_filter_loss_db = 1.0
if_amp_gain_db = 23.0
if_amp_nf_db = 0.0  # Negligible
output_filter_loss_db = 1.0

# System
bandwidth_khz = 100.0
temperature_k = 300.0

# --- Step 1: Calculate Effective Isotropic Radiated Power (EIRP) ---
total_tx_loss_db = tx_antenna_loss_db + tx_filter_loss_db + tx_cable_loss_db
eirp_dbm = tx_power_dbm + tx_antenna_gain_db - total_tx_loss_db

# --- Step 2: Calculate Free Space Path Loss (FSPL) ---
# FSPL (dB) = 20*log10(d_km) + 20*log10(f_GHz) + 92.45
fspl_db = 20 * math.log10(distance_km) + 20 * math.log10(frequency_ghz) + 92.45

# --- Step 3: Calculate Received Signal Power (S) ---
# This is the signal power at the receiver's input terminals (after the Rx antenna)
received_signal_power_dbm = eirp_dbm - fspl_db + rx_antenna_gain_db

# --- Step 4: Calculate Receiver's Cascaded Noise Figure (NF_rx) ---
# The receiver chain starts with the antenna loss.
# Format: (Gain in dB, Noise Figure in dB)
# For passive components (losses), Gain = -Loss and NF = Loss.
rx_stages = [
    (-rx_antenna_loss_db, rx_antenna_loss_db),
    (-rx_input_filter_loss_db, rx_input_filter_loss_db),
    (lna_gain_db, lna_nf_db),
    (-mixer_loss_db, mixer_loss_db),
    (-if_filter_loss_db, if_filter_loss_db),
    (if_amp_gain_db, if_amp_nf_db),
    (-output_filter_loss_db, output_filter_loss_db)
]

total_noise_factor = 0
cascaded_gain_linear = 1.0

for gain_db, nf_db in rx_stages:
    gain_linear = db_to_linear(gain_db)
    noise_factor_linear = db_to_linear(nf_db)
    
    # Friis formula: F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
    stage_contribution = (noise_factor_linear - 1) / cascaded_gain_linear
    total_noise_factor += stage_contribution
    
    cascaded_gain_linear *= gain_linear

# The first stage's noise factor is added differently (F1, not (F1-1)/G0)
# Our loop calculates (F-1)/G, so we need to adjust.
# F_total = (F1-1)/1 + (F2-1)/G1 + ... = (sum of terms)
# The actual total noise factor is F_total = F1 + (F2-1)/G1 + ...
# So, F_total = (sum of terms) + 1
total_noise_factor += 1 # This accounts for the first stage correctly
rx_nf_db = linear_to_db(total_noise_factor)

# --- Step 5: Calculate Total Input-Referred Noise Power (N) ---
# Thermal noise floor in dBm = -174 dBm/Hz + 10*log10(Bandwidth in Hz)
bandwidth_hz = bandwidth_khz * 1000
thermal_noise_floor_dbm = -174 + 10 * math.log10(bandwidth_hz)

# Total noise is thermal noise plus the noise added by the receiver
total_input_noise_dbm = thermal_noise_floor_dbm + rx_nf_db

# --- Step 6: Calculate Final SNR ---
snr_db = received_signal_power_dbm - total_input_noise_dbm

# --- Print the results ---
print("--- Link Budget Calculation ---")
print("\n1. Signal Power Calculation:")
print(f"  EIRP = Tx Power ({tx_power_dbm:.2f} dBm) + Tx Gain ({tx_antenna_gain_db:.2f} dB) - Tx Loss ({total_tx_loss_db:.2f} dB) = {eirp_dbm:.2f} dBm")
print(f"  Free Space Path Loss (FSPL) = {fspl_db:.2f} dB")
print(f"  Received Signal Power (S) = EIRP ({eirp_dbm:.2f} dBm) - FSPL ({fspl_db:.2f} dB) + Rx Gain ({rx_antenna_gain_db:.2f} dB)")
print(f"  S = {received_signal_power_dbm:.2f} dBm")

print("\n2. Noise Power Calculation:")
print(f"  Receiver Cascaded Noise Figure (NF_rx) = {rx_nf_db:.2f} dB")
print(f"  Thermal Noise Floor (kTB) = {thermal_noise_floor_dbm:.2f} dBm (for {bandwidth_khz:.0f} kHz BW)")
print(f"  Total Input Noise Power (N) = Thermal Noise ({thermal_noise_floor_dbm:.2f} dBm) + NF_rx ({rx_nf_db:.2f} dB)")
print(f"  N = {total_input_noise_dbm:.2f} dBm")

print("\n3. Final SNR Calculation:")
print(f"  SNR (dB) = Signal Power (S) - Noise Power (N)")
print(f"  SNR (dB) = {received_signal_power_dbm:.2f} dBm - ({total_input_noise_dbm:.2f} dBm)")
print(f"  Resulting SNR = {snr_db:.2f} dB")

# The final answer in the required format
final_answer = snr_db
print(f"\n<<<{final_answer:.2f}>>>")