import math

# Define a function to convert dB values to linear scale for gain and noise factor
def db_to_linear(db_val):
    return 10**(db_val / 10.0)

# --- Given Parameters ---
# Transmitter side
tx_power_dbm = 30.0
tx_ant_gain_db = 20.0
tx_ant_loss_db = 1.0
tx_filter_loss_db = 1.0
tx_cable_loss_db = 1.0

# Path
freq_ghz = 24.0
distance_km = 10.0
signal_bw_hz = 100 * 1000 # 100 kHz in Hz

# Receiver side
rx_ant_gain_db = 1.0
rx_ant_loss_db = 0.5
rx_filter_loss_db = 1.0 # Input filter
lna_gain_db = 36.0
lna_nf_db = 2.0
mixer_loss_db = 9.0 # This is conversion loss, treated as both gain and NF
if_filter_loss_db = 1.0
# The IF amp and subsequent filter have negligible impact on NF due to high preceding gain

# Constants
k = 1.38e-23  # Boltzmann's constant in J/K
T = 300.0     # Ambient temperature in Kelvin

print("Step-by-step SNR Calculation:\n")

# 1. Calculate Transmitter EIRP
print("--- 1. Transmitter EIRP Calculation ---")
total_tx_loss_db = tx_ant_loss_db + tx_filter_loss_db + tx_cable_loss_db
eirp_dbm = tx_power_dbm + tx_ant_gain_db - total_tx_loss_db
print(f"EIRP = Tx Power ({tx_power_dbm} dBm) + Tx Antenna Gain ({tx_ant_gain_db} dB) - Total Tx Loss ({total_tx_loss_db:.1f} dB)")
print(f"EIRP = {eirp_dbm:.2f} dBm\n")


# 2. Calculate Free Space Path Loss (FSPL)
print("--- 2. Free Space Path Loss (FSPL) Calculation ---")
fspl_db = 20 * math.log10(distance_km) + 20 * math.log10(freq_ghz) + 92.45
print(f"FSPL = 20*log10({distance_km} km) + 20*log10({freq_ghz} GHz) + 92.45")
print(f"FSPL = {fspl_db:.2f} dB\n")


# 3. Calculate Received Signal Power (at Rx Antenna Terminals)
print("--- 3. Received Signal Power Calculation ---")
# Reference point: after Rx antenna gain, before Rx losses.
received_signal_power_dbm = eirp_dbm - fspl_db + rx_ant_gain_db
print(f"Signal Power = EIRP ({eirp_dbm:.2f} dBm) - FSPL ({fspl_db:.2f} dB) + Rx Antenna Gain ({rx_ant_gain_db:.1f} dB)")
print(f"Received Signal Power = {received_signal_power_dbm:.2f} dBm\n")


# 4. Calculate Total Receiver Noise Figure
print("--- 4. Receiver Noise Figure (NF) Calculation ---")
# The "System" starts at our reference point.
# First stage: Rx losses (Antenna Loss + Input Filter Loss)
stage1_loss_db = rx_ant_loss_db + rx_filter_loss_db
F1 = db_to_linear(stage1_loss_db) # Noise Factor of a passive loss is the loss itself
G1 = db_to_linear(-stage1_loss_db)

# Second stage: LNA
F2 = db_to_linear(lna_nf_db)
G2 = db_to_linear(lna_gain_db)

# Third stage: Mixer
F3 = db_to_linear(mixer_loss_db) # NF of passive mixer is its loss
G3 = db_to_linear(-mixer_loss_db)

# Fourth stage: IF Filter
F4 = db_to_linear(if_filter_loss_db)

# Calculate cascaded Noise Figure using Friis formula
# F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
F_total = F1 + (F2 - 1) / G1 + (F3 - 1) / (G1 * G2) + (F4 - 1) / (G1 * G2 * G3)
total_nf_db = 10 * math.log10(F_total)
print(f"Total System NF is calculated using the Friis formula for cascaded components:")
print(f"Receiver Losses ({stage1_loss_db:.1f} dB), LNA (NF {lna_nf_db:.1f} dB), Mixer (NF {mixer_loss_db:.1f} dB), etc.")
print(f"Total Receiver Noise Figure (NF) = {total_nf_db:.2f} dB\n")

# 5. Calculate Total Receiver Noise Power
print("--- 5. Receiver Noise Power Calculation ---")
noise_floor_dbm_per_hz = 10 * math.log10(k * T * 1000) # x1000 to convert W to mW for dBm
noise_bw_factor_db = 10 * math.log10(signal_bw_hz)
total_noise_power_dbm = noise_floor_dbm_per_hz + noise_bw_factor_db + total_nf_db
print(f"Noise Power = Thermal Noise Floor + 10*log10(Bandwidth) + Total NF")
print(f"Noise Power = ({noise_floor_dbm_per_hz:.2f} dBm/Hz) + ({noise_bw_factor_db:.2f} dB) + ({total_nf_db:.2f} dB)")
print(f"Total Noise Power = {total_noise_power_dbm:.2f} dBm\n")


# 6. Calculate Final SNR
print("--- 6. Final Signal-to-Noise Ratio (SNR) Calculation ---")
snr_db = received_signal_power_dbm - total_noise_power_dbm
print(f"SNR (dB) = Received Signal Power (dBm) - Total Noise Power (dBm)")
print(f"SNR = {received_signal_power_dbm:.2f} dBm - ({total_noise_power_dbm:.2f} dBm)")
print(f"Final SNR = {snr_db:.2f} dB")
<<<28.27>>>