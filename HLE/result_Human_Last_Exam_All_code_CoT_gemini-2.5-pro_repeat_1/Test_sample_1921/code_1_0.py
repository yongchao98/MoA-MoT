import math

# --- Step-by-step calculation for Microwave Link SNR ---

# --- Given Parameters ---
tx_power_dbm = 30.0
tx_ant_gain_db = 20.0
tx_ant_loss_db = 1.0
tx_filter_loss_db = 1.0
tx_cable_loss_db = 1.0
freq_ghz = 24.0
bw_khz = 100.0
distance_km = 10.0
rx_ant_gain_db = 1.0
rx_ant_loss_db = 0.5
rx_filter_loss_db = 1.0
lna_gain_db = 36.0
lna_nf_db = 2.0
mixer_loss_db = 9.0
temp_k = 300.0
k_boltzmann = 1.38e-23  # Boltzmann's constant

# --- 1. Calculate Transmitter EIRP (Effective Isotropic Radiated Power) ---
print("--- Step 1: Calculate Transmitter EIRP ---")
tx_total_loss_db = tx_ant_loss_db + tx_filter_loss_db + tx_cable_loss_db
eirp_dbm = tx_power_dbm + tx_ant_gain_db - tx_total_loss_db
print(f"Transmitter Losses = {tx_ant_loss_db} dB (ant) + {tx_filter_loss_db} dB (filter) + {tx_cable_loss_db} dB (cable) = {tx_total_loss_db:.1f} dB")
print(f"EIRP = {tx_power_dbm} dBm (Tx Power) + {tx_ant_gain_db} dB (Tx Gain) - {tx_total_loss_db:.1f} dB (Tx Losses) = {eirp_dbm:.1f} dBm\n")

# --- 2. Calculate Free Space Path Loss (FSPL) ---
print("--- Step 2: Calculate Free Space Path Loss (FSPL) ---")
# FSPL (dB) = 20*log10(distance_km) + 20*log10(freq_ghz) + 92.45
fspl_db = 20 * math.log10(distance_km) + 20 * math.log10(freq_ghz) + 92.45
print(f"FSPL = 20*log10({distance_km}) + 20*log10({freq_ghz}) + 92.45 = {fspl_db:.2f} dB\n")

# --- 3. Calculate Received Signal Power (at Rx antenna output) ---
print("--- Step 3: Calculate Received Signal Power ---")
# This is the power available at the terminals of the Rx antenna, before any receiver losses.
rx_signal_power_dbm = eirp_dbm - fspl_db + rx_ant_gain_db
print(f"Received Signal Power = {eirp_dbm:.1f} dBm (EIRP) - {fspl_db:.2f} dB (FSPL) + {rx_ant_gain_db} dB (Rx Gain) = {rx_signal_power_dbm:.2f} dBm\n")

# --- 4. Calculate Receiver System Noise Figure (NF) ---
print("--- Step 4: Calculate Receiver System Noise Figure ---")
# The system noise figure is the sum of losses before the LNA, plus the LNA's noise figure.
# The noise from components after the high-gain LNA is negligible.
# Let's verify the noise contribution from the mixer is small.
lna_nf_lin = 10**(lna_nf_db / 10)
lna_gain_lin = 10**(lna_gain_db / 10)
mixer_nf_lin = 10**(mixer_loss_db / 10) # NF of a passive mixer is its conversion loss
# NF of LNA+Mixer combo using Friis formula: F_total = F1 + (F2-1)/G1
lna_mixer_nf_lin = lna_nf_lin + (mixer_nf_lin - 1) / lna_gain_lin
lna_mixer_nf_db = 10 * math.log10(lna_mixer_nf_lin)
# The NF of the active chain is dominated by the LNA's NF.
active_chain_nf_db = lna_nf_db # We use the LNA's NF as a very close approximation (2.007 dB calculated vs 2.0 dB given)
pre_lna_losses_db = rx_ant_loss_db + rx_filter_loss_db
system_nf_db = pre_lna_losses_db + active_chain_nf_db
print(f"Pre-LNA Losses = {rx_ant_loss_db} dB (Rx Ant) + {rx_filter_loss_db} dB (Rx Filter) = {pre_lna_losses_db:.1f} dB")
print(f"Active Chain NF is dominated by the LNA's NF of {active_chain_nf_db:.1f} dB.")
print(f"Total System NF = {pre_lna_losses_db:.1f} dB (Pre-LNA Losses) + {active_chain_nf_db:.1f} dB (LNA NF) = {system_nf_db:.1f} dB\n")

# --- 5. Calculate Total Input-Referred Noise Power ---
print("--- Step 5: Calculate Total Input-Referred Noise Power ---")
bw_hz = bw_khz * 1000
# Thermal Noise (N) = k * T * B
# In dBm: N_dBm = 10*log10(k*T*B*1000)
thermal_noise_power_dbm = 10 * math.log10(k_boltzmann * temp_k * bw_hz * 1000)
total_noise_power_dbm = thermal_noise_power_dbm + system_nf_db
print(f"Thermal Noise Power = 10*log10(k*T*B*1000) for B={bw_khz:.0f}kHz = {thermal_noise_power_dbm:.2f} dBm")
print(f"Total Noise Power = {thermal_noise_power_dbm:.2f} dBm (Thermal) + {system_nf_db:.1f} dB (System NF) = {total_noise_power_dbm:.2f} dBm\n")

# --- 6. Calculate Final Signal-to-Noise Ratio (SNR) ---
print("--- Step 6: Calculate Final SNR ---")
snr_db = rx_signal_power_dbm - total_noise_power_dbm
print(f"SNR (dB) = Received Signal ({rx_signal_power_dbm:.2f} dBm) - Total Noise ({total_noise_power_dbm:.2f} dBm)")
print(f"Final SNR = {rx_signal_power_dbm:.2f} - ({total_noise_power_dbm:.2f}) = {snr_db:.2f} dB")

# --- Final Answer ---
final_snr = round(snr_db, 2)