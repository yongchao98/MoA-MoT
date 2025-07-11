import math

def to_linear(db_value):
    """Converts a dB value to a linear scale."""
    return 10**(db_value / 10)

def to_db(linear_value):
    """Converts a linear value to dB."""
    return 10 * math.log10(linear_value)

# --- Given Parameters ---
# Transmitter (Tx)
tx_power_dBm = 30.0
tx_ant_gain_dB = 20.0
tx_ant_loss_dB = 1.0
tx_filter_loss_dB = 1.0
tx_cable_loss_dB = 1.0

# Path
frequency_GHz = 24.0
distance_km = 10.0

# Receiver (Rx)
rx_ant_gain_dB = 1.0
rx_ant_loss_dB = 0.5
rx_input_filter_loss_dB = 1.0
lna_gain_dB = 36.0
lna_nf_dB = 2.0
mixer_loss_dB = 9.0
mixer_nf_dB = 9.0  # Assuming Noise Figure equals conversion loss for a passive mixer
if_filter_loss_dB = 1.0
# The IF amplifier has negligible NF and the final output filter comes after,
# so they don't impact the SNR calculation.

# System
signal_bw_Hz = 100_000.0  # 100 kHz
temp_K = 300.0
k_boltzmann = 1.38e-23  # Boltzmann's constant in J/K

# --- Step 1: Calculate Transmitter EIRP ---
tx_total_loss_dB = tx_ant_loss_dB + tx_filter_loss_dB + tx_cable_loss_dB
eirp_dBm = tx_power_dBm + tx_ant_gain_dB - tx_total_loss_dB
print("Step 1: Calculate Transmitter EIRP")
print(f"EIRP (dBm) = Tx Power (dBm) + Tx Antenna Gain (dB) - Tx Total Loss (dB)")
print(f"EIRP (dBm) = {tx_power_dBm} + {tx_ant_gain_dB} - ({tx_ant_loss_dB} + {tx_filter_loss_dB} + {tx_cable_loss_dB})")
print(f"EIRP (dBm) = {eirp_dBm:.2f}\n")

# --- Step 2: Calculate Free Space Path Loss (FSPL) ---
fspl_dB = 20 * math.log10(distance_km) + 20 * math.log10(frequency_GHz) + 92.45
print("Step 2: Calculate Free Space Path Loss (FSPL)")
print(f"FSPL (dB) = 20*log10(distance_km) + 20*log10(frequency_GHz) + 92.45")
print(f"FSPL (dB) = 20*log10({distance_km}) + 20*log10({frequency_GHz}) + 92.45")
print(f"FSPL (dB) = {fspl_dB:.2f}\n")

# --- Step 3: Calculate Received Signal Power (Prx) ---
prx_dBm = eirp_dBm - fspl_dB + rx_ant_gain_dB
print("Step 3: Calculate Received Signal Power (Prx)")
print(f"Prx (dBm) = EIRP (dBm) - FSPL (dB) + Rx Antenna Gain (dB)")
print(f"Prx (dBm) = {eirp_dBm:.2f} - {fspl_dB:.2f} + {rx_ant_gain_dB}")
print(f"Prx (dBm) = {prx_dBm:.2f}\n")

# --- Step 4: Calculate Total Receiver Noise Figure ---
# Convert gains and NFs of the main chain (post-losses) to linear scale
lna_gain_lin = to_linear(lna_gain_dB)
lna_nf_lin = to_linear(lna_nf_dB)
mixer_gain_lin = to_linear(-mixer_loss_dB) # Loss is negative gain
mixer_nf_lin = to_linear(mixer_nf_dB)
if_filter_gain_lin = to_linear(-if_filter_loss_dB)
if_filter_nf_lin = to_linear(if_filter_loss_dB) # NF of a passive filter equals its loss

# Use Friis formula for the cascaded noise figure of the LNA, Mixer, and IF Filter
# F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
nf_chain_lin = lna_nf_lin + (mixer_nf_lin - 1) / lna_gain_lin + (if_filter_nf_lin - 1) / (lna_gain_lin * mixer_gain_lin)
nf_chain_dB = to_db(nf_chain_lin)

# Add losses before the LNA (these add directly to the NF in dB)
pre_lna_loss_dB = rx_ant_loss_dB + rx_input_filter_loss_dB
total_nf_dB = pre_lna_loss_dB + nf_chain_dB
print("Step 4: Calculate Total Receiver Noise Figure (NF)")
print(f"NF_chain (dB) = NF of LNA, Mixer, and IF Filter = {nf_chain_dB:.2f} dB")
print(f"Total NF (dB) = Pre-LNA Losses (dB) + NF_chain (dB)")
print(f"Total NF (dB) = ({rx_ant_loss_dB} + {rx_input_filter_loss_dB}) + {nf_chain_dB:.2f}")
print(f"Total NF (dB) = {total_nf_dB:.2f}\n")

# --- Step 5: Calculate Input-Referred Noise Power ---
# Thermal noise in Watts: N = k * T * B
# Thermal noise in dBm: N = 10*log10(k * T * B * 1000)
noise_thermal_dBm = 10 * math.log10(k_boltzmann * temp_K * signal_bw_Hz * 1000)
noise_input_total_dBm = noise_thermal_dBm + total_nf_dB
print("Step 5: Calculate Input-Referred Noise Power")
print(f"Thermal Noise (dBm) = 10*log10(k * T * B * 1000)")
print(f"Thermal Noise (dBm) = 10*log10({k_boltzmann:.2e} * {temp_K} * {signal_bw_Hz} * 1000) = {noise_thermal_dBm:.2f} dBm")
print(f"Total Input Noise (dBm) = Thermal Noise (dBm) + Total NF (dB)")
print(f"Total Input Noise (dBm) = {noise_thermal_dBm:.2f} + {total_nf_dB:.2f}")
print(f"Total Input Noise (dBm) = {noise_input_total_dBm:.2f}\n")

# --- Step 6: Calculate Final SNR ---
snr_dB = prx_dBm - noise_input_total_dBm
print("Step 6: Calculate Final SNR")
print(f"SNR (dB) = Received Signal Power (dBm) - Total Input Noise (dBm)")
print(f"SNR (dB) = {prx_dBm:.2f} - ({noise_input_total_dBm:.2f})")
print(f"SNR (dB) = {snr_dB:.2f}")

final_snr = round(snr_dB, 2)