import math

# This script calculates the Signal-to-Noise Ratio (SNR) for a microwave link.

# --- Helper Functions ---
def db_to_lin(db_val):
    """Converts a value from dB to linear scale (for gain/loss)."""
    return 10**(db_val / 10)

def lin_to_db(lin_val):
    """Converts a value from linear scale to dB."""
    return 10 * math.log10(lin_val)

def lin_to_dbm(lin_val_watts):
    """Converts a power value from Watts to dBm."""
    return 10 * math.log10(lin_val_watts / 0.001)

# --- Constants ---
k = 1.38e-23  # Boltzmann's constant (J/K)
T = 300       # Ambient temperature (K)

# --- System Parameters ---
# Transmitter (Tx)
tx_power_dbm = 30.0
tx_ant_gain_db = 20.0
tx_ant_loss_db = 1.0
tx_filter_loss_db = 1.0
tx_cable_loss_db = 1.0
frequency_ghz = 24.0
bandwidth_khz = 100.0

# Path
distance_km = 10.0

# Receiver (Rx)
rx_ant_gain_db = 1.0
rx_ant_loss_db = 0.5
rx_filter_loss_db = 1.0
lna_gain_db = 36.0
lna_nf_db = 2.0
mixer_loss_db = 9.0
if_filter_loss_db = 1.0
if_amp_gain_db = 23.0
if_amp_nf_db = 0.0 # Negligible
output_filter_loss_db = 1.0

# --- Calculations ---

# 1. Effective Isotropic Radiated Power (EIRP)
print("--- 1. Transmitter EIRP Calculation ---")
tx_total_loss_db = tx_ant_loss_db + tx_filter_loss_db + tx_cable_loss_db
eirp_dbm = tx_power_dbm + tx_ant_gain_db - tx_total_loss_db
print(f"EIRP (dBm) = Tx Power (dBm) + Tx Ant Gain (dB) - Tx Losses (dB)")
print(f"EIRP (dBm) = {tx_power_dbm} + {tx_ant_gain_db} - ({tx_ant_loss_db} + {tx_filter_loss_db} + {tx_cable_loss_db}) = {eirp_dbm:.2f} dBm\n")

# 2. Free Space Path Loss (FSPL)
print("--- 2. Free Space Path Loss (FSPL) Calculation ---")
fspl_db = 20 * math.log10(distance_km) + 20 * math.log10(frequency_ghz) + 92.45
print(f"FSPL (dB) = 20*log10(d_km) + 20*log10(f_GHz) + 92.45")
print(f"FSPL (dB) = 20*log10({distance_km}) + 20*log10({frequency_ghz}) + 92.45 = {fspl_db:.2f} dB\n")

# 3. Received Signal Power at LNA Input (P_signal)
print("--- 3. Received Signal Power Calculation ---")
rx_frontend_losses_db = rx_ant_loss_db + rx_filter_loss_db
signal_power_at_lna_input_dbm = eirp_dbm - fspl_db + rx_ant_gain_db - rx_frontend_losses_db
print(f"P_signal (dBm) = EIRP (dBm) - FSPL (dB) + Rx Ant Gain (dB) - Rx Frontend Losses (dB)")
print(f"P_signal (dBm) = {eirp_dbm:.2f} - {fspl_db:.2f} + {rx_ant_gain_db} - ({rx_ant_loss_db} + {rx_filter_loss_db}) = {signal_power_at_lna_input_dbm:.2f} dBm\n")

# 4. Receiver Noise Figure (NF)
print("--- 4. Receiver Noise Calculation ---")
# Gains in linear scale
g1_lna = db_to_lin(lna_gain_db)
g2_mixer = db_to_lin(-mixer_loss_db) # Loss is negative gain
g3_if_filter = db_to_lin(-if_filter_loss_db)
g4_if_amp = db_to_lin(if_amp_gain_db)

# Noise Factors (linear scale)
f1_lna = db_to_lin(lna_nf_db)
f2_mixer = db_to_lin(mixer_loss_db) # For passive components, NF = Loss
f3_if_filter = db_to_lin(if_filter_loss_db)
f4_if_amp = db_to_lin(if_amp_nf_db) # Negligible NF = 0 dB -> F = 1
f5_output_filter = db_to_lin(output_filter_loss_db)

# Friis formula for cascaded noise factor: F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
f_chain_lin = f1_lna + (f2_mixer - 1) / g1_lna + \
              (f3_if_filter - 1) / (g1_lna * g2_mixer) + \
              (f4_if_amp - 1) / (g1_lna * g2_mixer * g3_if_filter) + \
              (f5_output_filter - 1) / (g1_lna * g2_mixer * g3_if_filter * g4_if_amp)

nf_chain_db = lin_to_db(f_chain_lin)
print(f"Cascaded Noise Figure of Rx Chain (from LNA) = {nf_chain_db:.2f} dB")

# Total system noise figure includes frontend losses
nf_system_db = rx_frontend_losses_db + nf_chain_db
print(f"Total System Noise Figure (dB) = Frontend Losses (dB) + Rx Chain NF (dB)")
print(f"Total System Noise Figure (dB) = {rx_frontend_losses_db} + {nf_chain_db:.2f} = {nf_system_db:.2f} dB\n")

# 5. Total Noise Power (P_noise)
print("--- 5. Total Noise Power Calculation ---")
bandwidth_hz = bandwidth_khz * 1000
thermal_noise_watts = k * T * bandwidth_hz
thermal_noise_dbm = lin_to_dbm(thermal_noise_watts)

total_noise_power_dbm = thermal_noise_dbm + nf_system_db
print(f"P_noise (dBm) = Thermal Noise (kTB) (dBm) + System NF (dB)")
print(f"P_noise (dBm) = {thermal_noise_dbm:.2f} + {nf_system_db:.2f} = {total_noise_power_dbm:.2f} dBm\n")

# 6. Signal-to-Noise Ratio (SNR)
print("--- 6. Final SNR Calculation ---")
snr_db = signal_power_at_lna_input_dbm - total_noise_power_dbm
print(f"SNR (dB) = P_signal (dBm) - P_noise (dBm)")
print(f"SNR (dB) = {signal_power_at_lna_input_dbm:.2f} - ({total_noise_power_dbm:.2f}) = {snr_db:.2f} dB")

# Final Answer
print("\nThe final resulting SNR of the received signal is:")
print(f"{snr_db:.2f} dB")
<<<26.77>>>