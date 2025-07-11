import math

# A function to convert dB values to linear scale for calculations
def db_to_linear(db_val):
    """Converts a value in dB to a linear scale."""
    return 10**(db_val / 10)

# A function to convert linear values back to dB
def linear_to_db(lin_val):
    """Converts a linear value to dB scale."""
    return 10 * math.log10(lin_val)

# --- Input Parameters from the problem description ---

# Transmitter Parameters
tx_power_dBm = 30.0
tx_ant_gain_dB = 20.0
tx_ant_loss_dB = 1.0
tx_filter_loss_dB = 1.0
tx_cable_loss_dB = 1.0

# Path Parameters
frequency_GHz = 24.0
distance_km = 10.0

# Receiver Parameters (in order of signal path)
rx_ant_gain_dB = 1.0
rx_ant_loss_dB = 0.5
rx_input_filter_loss_dB = 1.0
lna_gain_dB = 36.0
lna_nf_dB = 2.0
mixer_loss_dB = 9.0
mixer_nf_dB = 9.0  # For a passive component like a mixer, NF is equal to its loss
if_filter_loss_dB = 1.0
if_amp_gain_dB = 23.0
if_amp_nf_dB = 0.0 # Negligible NF is treated as 0 dB
output_filter_loss_dB = 1.0

# System Parameters
signal_bw_Hz = 100000.0  # 100 kHz

print("### Calculating Signal-to-Noise Ratio (SNR) ###\n")

# Step 1: Calculate Effective Isotropic Radiated Power (EIRP)
print("--- Step 1: Calculate Signal Power ---")
eirp_dBm = tx_power_dBm + tx_ant_gain_dB - tx_ant_loss_dB - tx_filter_loss_dB - tx_cable_loss_dB
print(f"Transmitter EIRP = {tx_power_dBm} dBm (Tx) + {tx_ant_gain_dB} dB (Ant Gain) - {tx_ant_loss_dB} dB (Ant Loss) - {tx_filter_loss_dB} dB (Filter) - {tx_cable_loss_dB} dB (Cable)")
print(f"Calculated EIRP = {eirp_dBm:.2f} dBm\n")

# Step 2: Calculate Free Space Path Loss (FSPL)
fspl_dB = 20 * math.log10(distance_km) + 20 * math.log10(frequency_GHz) + 92.45
print(f"Free Space Path Loss (FSPL) for {distance_km} km at {frequency_GHz} GHz = {fspl_dB:.2f} dB\n")

# Step 3: Calculate Received Signal Power at the LNA input
# Losses before the LNA include the receiver antenna loss and the input filter loss.
rx_pre_lna_losses_dB = rx_ant_loss_dB + rx_input_filter_loss_dB
rx_signal_power_dBm = eirp_dBm - fspl_dB + rx_ant_gain_dB - rx_pre_lna_losses_dB
print("Received Signal Power = EIRP - FSPL + Rx Ant Gain - Rx Losses (before LNA)")
print(f"Received Signal Power = {eirp_dBm:.2f} dBm - {fspl_dB:.2f} dB + {rx_ant_gain_dB} dB - {rx_pre_lna_losses_dB} dB")
print(f"Calculated Received Signal Power = {rx_signal_power_dBm:.2f} dBm\n")


# Step 4: Calculate Total Noise Power at the LNA input
print("--- Step 2: Calculate Noise Power ---")
# Thermal Noise Power using the standard formula for 300K temperature (-174 dBm/Hz)
thermal_noise_power_dBm = -174 + 10 * math.log10(signal_bw_Hz)
print(f"Thermal Noise in {signal_bw_Hz / 1000} kHz Bandwidth = -174 dBm/Hz + 10*log10({signal_bw_Hz}) Hz = {thermal_noise_power_dBm:.2f} dBm\n")

# Calculate Receiver System Noise Figure using Friis Formula
# We calculate it from the LNA onwards, as our signal/noise reference point is the LNA input.
# Convert necessary dB values to linear scale for the formula
g_lna_lin = db_to_linear(lna_gain_dB)
f_lna_lin = db_to_linear(lna_nf_dB)

g_mixer_lin = db_to_linear(-mixer_loss_dB) # Gain is negative for a loss
f_mixer_lin = db_to_linear(mixer_nf_dB)

g_if_filt_lin = db_to_linear(-if_filter_loss_dB)
f_if_filt_lin = db_to_linear(if_filter_loss_dB) # NF of passive filter is its loss

# Friis formula: F_sys = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
stage1_f = f_lna_lin
stage2_contrib = (f_mixer_lin - 1) / g_lna_lin
stage3_contrib = (f_if_filt_lin - 1) / (g_lna_lin * g_mixer_lin)

f_sys_lin = stage1_f + stage2_contrib + stage3_contrib
# The remaining stages have negligible contribution due to the high gain of the LNA
# and the division by multiple gain/loss stages.

nf_sys_dB = linear_to_db(f_sys_lin)
print("Receiver System Noise Figure (from LNA onwards):")
print(f"  Noise Factor = F_lna + (F_mixer-1)/G_lna + (F_if_filter-1)/(G_lna*G_mixer) + ...")
print(f"  Noise Factor = {stage1_f:.3f} + {stage2_contrib:.3f} + {stage3_contrib:.3f} = {f_sys_lin:.3f}")
print(f"Calculated System Noise Figure = {nf_sys_dB:.2f} dB\n")

# Total Noise Power = Thermal Noise + System Noise Figure
total_noise_power_dBm = thermal_noise_power_dBm + nf_sys_dB
print("Total Noise Power (at LNA input) = Thermal Noise Power + System Noise Figure")
print(f"Total Noise Power = {thermal_noise_power_dBm:.2f} dBm + {nf_sys_dB:.2f} dB")
print(f"Calculated Total Noise Power = {total_noise_power_dBm:.2f} dBm\n")

# Step 5: Calculate Final SNR
print("--- Step 3: Calculate Final SNR ---")
snr_dB = rx_signal_power_dBm - total_noise_power_dBm
print("SNR (dB) = Received Signal Power (dBm) - Total Noise Power (dBm)")
print(f"Final SNR (dB) = {rx_signal_power_dBm:.2f} - ({total_noise_power_dBm:.2f})")
print(f"Resulting SNR = {snr_dB:.2f} dB")

# Final answer in the required format
final_answer = round(snr_dB, 2)
print(f"\n<<<{final_answer}>>>")
