import math

# This script calculates the Signal-to-Noise Ratio (SNR) for a microwave link.

# Plan:
# 1. Define all constants and parameters from the problem description.
# 2. Calculate the transmitter's Effective Isotropic Radiated Power (EIRP).
# 3. Calculate the Free Space Path Loss (FSPL).
# 4. Calculate the signal power (S) at the input of the Low Noise Amplifier (LNA).
# 5. Calculate the receiver's system noise figure (NF) referenced to the LNA input.
# 6. Calculate the thermal noise power floor.
# 7. Calculate the total noise power (N) at the LNA input.
# 8. Calculate the Signal-to-Noise Ratio (SNR) by subtracting N from S.
# 9. Print a breakdown of the calculation and the final result.

# Step 1: Define constants and parameters
# Constants
k = 1.380649e-23  # Boltzmann's constant (J/K)
T = 300.0         # Ambient temperature (K)
f = 24e9          # Frequency (Hz)
B = 100e3         # Signal Bandwidth (Hz)
d = 10e3          # Distance (m)
c = 299792458.0   # Speed of light (m/s)

# Transmitter (Tx) parameters in dB or dBm
tx_power_dBm = 30.0
tx_cable_loss_dB = 1.0
tx_filter_loss_dB = 1.0
tx_antenna_loss_dB = 1.0
tx_antenna_gain_dBi = 20.0

# Receiver (Rx) parameters in dB
rx_antenna_gain_dBi = 1.0
rx_antenna_loss_dB = 0.5
rx_filter_loss_dB = 1.0
lna_gain_dB = 36.0
lna_nf_dB = 2.0
mixer_loss_dB = 9.0      # This is conversion loss, Gain = -9 dB
mixer_nf_dB = 9.0        # For a passive mixer, NF = Loss
if_filter_loss_dB = 1.0
if_amp_nf_dB = 0.0       # Negligible NF
if_amp_gain_dB = 23.0
out_filter_loss_dB = 1.0

# Helper functions for converting between dB and linear scale
def to_linear(db_value):
    """Converts a dB value to a linear ratio."""
    return 10**(db_value / 10.0)

def to_db(linear_value):
    """Converts a linear ratio to a dB value."""
    return 10 * math.log10(linear_value)

# Step 2: Calculate EIRP
tx_total_loss_dB = tx_cable_loss_dB + tx_filter_loss_dB + tx_antenna_loss_dB
eirp_dBm = tx_power_dBm - tx_total_loss_dB + tx_antenna_gain_dBi

# Step 3: Calculate Free Space Path Loss (FSPL)
fspl_dB = to_db((4 * math.pi * d * f / c)**2)

# Step 4: Calculate Signal Power (S) at the LNA input
losses_before_lna_dB = rx_antenna_loss_dB + rx_filter_loss_dB
signal_at_lna_input_dBm = eirp_dBm - fspl_dB + rx_antenna_gain_dBi - losses_before_lna_dB

# Step 5: Calculate cascaded Noise Figure (NF) at the LNA input using Friis formula
# Friis Formula for Noise Factor: F_sys = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
F_lna = to_linear(lna_nf_dB)
G_lna = to_linear(lna_gain_dB)
F_mixer = to_linear(mixer_nf_dB)
G_mixer = to_linear(-mixer_loss_dB) # Gain is the negative of loss
F_if_filter = to_linear(if_filter_loss_dB) # NF of a passive component is its loss
G_if_filter = to_linear(-if_filter_loss_dB)
F_if_amp = to_linear(if_amp_nf_dB)
G_if_amp = to_linear(if_amp_gain_dB)
F_out_filter = to_linear(out_filter_loss_dB)

# Calculate total noise factor referenced to LNA input
F_sys_linear = F_lna + \
               (F_mixer - 1) / G_lna + \
               (F_if_filter - 1) / (G_lna * G_mixer) + \
               (F_if_amp - 1) / (G_lna * G_mixer * G_if_filter) + \
               (F_out_filter - 1) / (G_lna * G_mixer * G_if_filter * G_if_amp)
nf_sys_at_lna_input_dB = to_db(F_sys_linear)

# Step 6: Calculate thermal noise power floor in dBm
thermal_noise_watts = k * T * B
thermal_noise_dBm = to_db(thermal_noise_watts) + 30

# Step 7: Calculate total Noise Power (N) at the LNA input
noise_at_lna_input_dBm = thermal_noise_dBm + nf_sys_at_lna_input_dB

# Step 8: Calculate the final SNR
snr_dB = signal_at_lna_input_dBm - noise_at_lna_input_dBm

# Step 9: Print the final equation and the result
print("--- Link Budget & SNR Calculation ---")

# Print Signal Calculation
print("\n--- Signal Power (S) Calculation ---")
print(f"EIRP = Tx Power ({tx_power_dBm:.1f} dBm) - Tx Losses ({tx_total_loss_dB:.1f} dB) + Tx Gain ({tx_antenna_gain_dBi:.1f} dBi) = {eirp_dBm:.2f} dBm")
print(f"Path Loss (FSPL) = {fspl_dB:.2f} dB")
print(f"Rx Losses before LNA = {losses_before_lna_dB:.1f} dB")
print(f"S at LNA input = EIRP ({eirp_dBm:.2f} dBm) - FSPL ({fspl_dB:.2f} dB) + Rx Gain ({rx_antenna_gain_dBi:.1f} dBi) - Rx Losses ({losses_before_lna_dB:.1f} dB) = {signal_at_lna_input_dBm:.2f} dBm")

# Print Noise Calculation
print("\n--- Noise Power (N) Calculation ---")
print(f"Thermal Noise in {B/1000:.0f} kHz BW = {thermal_noise_dBm:.2f} dBm")
print(f"Receiver System Noise Figure = {nf_sys_at_lna_input_dB:.2f} dB")
print(f"N at LNA input = Thermal Noise ({thermal_noise_dBm:.2f} dBm) + System NF ({nf_sys_at_lna_input_dB:.2f} dB) = {noise_at_lna_input_dBm:.2f} dBm")

# Print Final SNR Equation
print("\n--- Final SNR Calculation ---")
print(f"SNR (dB) = Signal Power (dBm) - Noise Power (dBm)")
print(f"SNR (dB) = {signal_at_lna_input_dBm:.2f} - ({noise_at_lna_input_dBm:.2f})")
print(f"\nResulting SNR = {snr_dB:.2f} dB")

print(f"\n<<<{snr_dB:.2f}>>>")