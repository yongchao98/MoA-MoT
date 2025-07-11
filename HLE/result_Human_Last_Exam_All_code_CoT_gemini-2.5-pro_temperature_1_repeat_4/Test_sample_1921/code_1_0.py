import math

# --- Given Parameters ---

# Transmitter (Tx)
tx_power_dbm = 30.0
tx_antenna_gain_db = 20.0
tx_antenna_loss_db = 1.0
tx_filter_loss_db = 1.0
tx_cable_loss_db = 1.0
freq_hz = 24e9  # 24 GHz
bandwidth_hz = 100e3 # 100 kHz

# Path
distance_m = 10e3 # 10 km

# Receiver (Rx)
rx_antenna_gain_db = 1.0
rx_antenna_loss_db = 0.5
rx_input_filter_loss_db = 1.0
lna_gain_db = 36.0
lna_nf_db = 2.0
mixer_loss_db = 9.0  # Conversion Loss
mixer_nf_db = 9.0    # Mixer NF is equal to its conversion loss
if_filter_loss_db = 1.0
if_amp_gain_db = 23.0
if_amp_nf_db = 0.0    # Negligible NF
output_filter_loss_db = 1.0

# Constants
k = 1.380649e-23  # Boltzmann's constant in J/K
T = 300.0         # Ambient temperature in Kelvin
c = 299792458.0   # Speed of light in m/s

# --- Calculations ---

# 1. Effective Isotropic Radiated Power (EIRP)
print("--- 1. Calculating Transmitter EIRP ---")
eirp_dbm = tx_power_dbm + tx_antenna_gain_db - tx_antenna_loss_db - tx_filter_loss_db - tx_cable_loss_db
print(f"EIRP (dBm) = Tx Power ({tx_power_dbm} dBm) + Tx Ant Gain ({tx_antenna_gain_db} dB) - Tx Losses ({tx_antenna_loss_db + tx_filter_loss_db + tx_cable_loss_db} dB)")
print(f"EIRP = {eirp_dbm:.2f} dBm\n")


# 2. Free Space Path Loss (FSPL)
print("--- 2. Calculating Free Space Path Loss (FSPL) ---")
# FSPL (dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
fspl_db = 20 * math.log10(distance_m) + 20 * math.log10(freq_hz) + 20 * math.log10(4 * math.pi / c)
print(f"FSPL (dB) = 20*log10({distance_m:.0f} m) + 20*log10({freq_hz:.0f} Hz) + 20*log10(4*pi/c)")
print(f"FSPL = {fspl_db:.2f} dB\n")


# 3. Received Signal Power (at Rx antenna port)
print("--- 3. Calculating Received Signal Power ---")
received_power_dbm = eirp_dbm - fspl_db + rx_antenna_gain_db
print(f"Received Power (dBm) = EIRP ({eirp_dbm:.2f} dBm) - FSPL ({fspl_db:.2f} dB) + Rx Ant Gain ({rx_antenna_gain_db} dB)")
print(f"Received Power = {received_power_dbm:.2f} dBm\n")


# 4. Total System Noise Figure (NF)
print("--- 4. Calculating System Noise Figure ---")
# Losses before the LNA add directly to the noise figure
loss_pre_lna_db = rx_antenna_loss_db + rx_input_filter_loss_db
print(f"Losses before LNA = Rx Ant Loss ({rx_antenna_loss_db} dB) + Rx Filter Loss ({rx_input_filter_loss_db} dB) = {loss_pre_lna_db:.2f} dB")

# Convert gains and NFs to linear scale for Friis formula
lna_gain = 10**(lna_gain_db / 10)
lna_nf_factor = 10**(lna_nf_db / 10)

mixer_gain = 10**(-mixer_loss_db / 10) # loss is negative gain
mixer_nf_factor = 10**(mixer_nf_db / 10)

if_filter_gain = 10**(-if_filter_loss_db / 10)
if_filter_nf_factor = 10**(if_filter_loss_db / 10) # For a passive device, NF factor = Loss factor

if_amp_gain = 10**(if_amp_gain_db / 10)
if_amp_nf_factor = 10**(if_amp_nf_db / 10)

# Friis formula for noise factor of the cascade starting from the LNA
nf_cascade_factor = (
    lna_nf_factor +
    (mixer_nf_factor - 1) / lna_gain +
    (if_filter_nf_factor - 1) / (lna_gain * mixer_gain) +
    (if_amp_nf_factor - 1) / (lna_gain * mixer_gain * if_filter_gain)
)
nf_cascade_db = 10 * math.log10(nf_cascade_factor)
print(f"Noise Figure of the active cascade (LNA onwards) = {nf_cascade_db:.2f} dB")

# Total system NF is the pre-LNA loss plus the cascade NF
nf_system_db = loss_pre_lna_db + nf_cascade_db
print(f"Total System NF (dB) = Pre-LNA Loss ({loss_pre_lna_db:.2f} dB) + Cascade NF ({nf_cascade_db:.2f} dB)")
print(f"Total System NF = {nf_system_db:.2f} dB\n")


# 5. Thermal Noise Floor
print("--- 5. Calculating Thermal Noise Floor ---")
thermal_noise_watts = k * T * bandwidth_hz
thermal_noise_dbm = 10 * math.log10(thermal_noise_watts * 1000) # Convert W to mW for dBm
print(f"Thermal Noise (dBm) = 10*log10(k * T * B * 1000) = 10*log10({k:.3e} * {T:.0f} * {bandwidth_hz:.0f} * 1000)")
print(f"Thermal Noise = {thermal_noise_dbm:.2f} dBm\n")


# 6. Total Input-Referred Noise
print("--- 6. Calculating Total Input-Referred Noise ---")
total_input_noise_dbm = thermal_noise_dbm + nf_system_db
print(f"Total Input Noise (dBm) = Thermal Noise ({thermal_noise_dbm:.2f} dBm) + System NF ({nf_system_db:.2f} dB)")
print(f"Total Input Noise = {total_input_noise_dbm:.2f} dBm\n")


# 7. Final Signal-to-Noise Ratio (SNR)
print("--- 7. Final SNR Calculation ---")
snr_db = received_power_dbm - total_input_noise_dbm
print("SNR (dB) = Received Signal Power (dBm) - Total Input Noise (dBm)")
print(f"SNR (dB) = {received_power_dbm:.2f} dBm - ({total_input_noise_dbm:.2f} dBm)")
print(f"SNR (dB) = {received_power_dbm:.2f} - {total_input_noise_dbm:.2f}")
print(f"Final SNR = {snr_db:.2f} dB")

<<<28.22>>>