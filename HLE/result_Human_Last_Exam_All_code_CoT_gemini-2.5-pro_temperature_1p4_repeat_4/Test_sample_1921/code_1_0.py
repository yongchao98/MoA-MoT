import math

# Step-by-step calculation of the SNR for the microwave link.

# --- 1. Transmitter (Tx) Side Calculations ---
tx_power_dbm = 30  # dBm
tx_antenna_gain_db = 20  # dB
tx_antenna_loss_db = 1  # dB
tx_filter_loss_db = 1  # dB
tx_cable_loss_db = 1  # dB

# Total Tx losses
total_tx_loss_db = tx_antenna_loss_db + tx_filter_loss_db + tx_cable_loss_db

# EIRP = Tx Power + Tx Antenna Gain - Total Tx Losses
eirp_dbm = tx_power_dbm + tx_antenna_gain_db - total_tx_loss_db
print(f"Calculating EIRP...")
print(f"EIRP = {tx_power_dbm} dBm + {tx_antenna_gain_db} dB - ({tx_antenna_loss_db} dB + {tx_filter_loss_db} dB + {tx_cable_loss_db} dB) = {eirp_dbm:.2f} dBm\n")


# --- 2. Path Loss Calculation ---
frequency_ghz = 24  # GHz
distance_km = 10  # km

# FSPL (dB) = 32.45 + 20*log10(f_MHz) + 20*log10(d_km)
frequency_mhz = frequency_ghz * 1000
fspl_db = 32.45 + 20 * math.log10(frequency_mhz) + 20 * math.log10(distance_km)
print(f"Calculating Free Space Path Loss (FSPL)...")
print(f"FSPL = 32.45 + 20*log10({frequency_mhz:.0f} MHz) + 20*log10({distance_km} km) = {fspl_db:.2f} dB\n")


# --- 3. Received Signal Power (S) Calculation ---
rx_antenna_gain_db = 1  # dB

# S = EIRP - FSPL + Rx Antenna Gain
signal_power_dbm = eirp_dbm - fspl_db + rx_antenna_gain_db
print(f"Calculating Received Signal Power (S)...")
print(f"S = {eirp_dbm:.2f} dBm - {fspl_db:.2f} dB + {rx_antenna_gain_db} dB = {signal_power_dbm:.2f} dBm\n")


# --- 4. Total Input-Referred Noise (N) Calculation ---
# --- 4a. Thermal Noise Floor ---
bandwidth_khz = 100  # kHz
bandwidth_hz = bandwidth_khz * 1000
thermal_noise_constant_dbm_hz = -174  # dBm/Hz at 300K

thermal_noise_dbm = thermal_noise_constant_dbm_hz + 10 * math.log10(bandwidth_hz)
print(f"Calculating Thermal Noise Floor...")
print(f"Thermal Noise = {thermal_noise_constant_dbm_hz} dBm/Hz + 10*log10({bandwidth_hz} Hz) = {thermal_noise_dbm:.2f} dBm\n")


# --- 4b. Receiver System Noise Figure (NF) ---
rx_antenna_loss_db = 0.5  # dB
rx_filter_loss_db = 1  # dB
lna_nf_db = 2  # dB

# System NF = Sum of losses before LNA + LNA NF
system_nf_db = rx_antenna_loss_db + rx_filter_loss_db + lna_nf_db
print(f"Calculating Receiver System Noise Figure (NF)...")
print(f"System NF = {rx_antenna_loss_db} dB (Ant Loss) + {rx_filter_loss_db} dB (Filter Loss) + {lna_nf_db} dB (LNA NF) = {system_nf_db:.2f} dB\n")


# --- 4c. Total Noise Power ---
total_noise_power_dbm = thermal_noise_dbm + system_nf_db
print(f"Calculating Total Input-Referred Noise Power (N)...")
print(f"N = {thermal_noise_dbm:.2f} dBm + {system_nf_db:.2f} dB = {total_noise_power_dbm:.2f} dBm\n")


# --- 5. Final SNR Calculation ---
snr_db = signal_power_dbm - total_noise_power_dbm
print(f"Final SNR Calculation:")
print(f"SNR (dB) = S (dBm) - N (dBm)")
# Final Answer printout
print(f"SNR (dB) = {signal_power_dbm:.2f} dBm - ({total_noise_power_dbm:.2f} dBm)")
print(f"Resulting SNR = {snr_db:.2f} dB")

# Return the final numerical answer in the required format
#<<<28.45>>>