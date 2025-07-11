import math

# --- Given Parameters ---
# Constants
k = 1.38e-23  # Boltzmann constant in J/K
T = 300       # Ambient temperature in Kelvin
c = 3e8       # Speed of light in m/s

# Transmitter (Tx)
P_tx_dBm = 30.0
G_tx_dB = 20.0
L_ant_tx_dB = 1.0
L_filt_tx_dB = 1.0
L_cable_tx_dB = 1.0
freq_Hz = 24e9
B_Hz = 100e3 # Signal Bandwidth

# Path
d_m = 10e3 # 10 km

# Receiver (Rx)
G_rx_dB = 1.0
L_ant_rx_dB = 0.5
L_filt1_rx_dB = 1.0
G_lna_dB = 36.0
NF_lna_dB = 2.0
L_mix_dB = 9.0
NF_mix_dB = 9.0 # Assume NF=Loss for passive mixer
L_filt_if_dB = 1.0
G_if_dB = 23.0
NF_if_dB = 0.0 # Negligible
L_filt_out_dB = 1.0

print("This script calculates the SNR for the given microwave link.")
print("The calculation is performed in dB/dBm and follows a standard link budget analysis.\n")

# --- Step 1: Calculate Effective Isotropic Radiated Power (EIRP) ---
print("--- 1. Signal Power Calculation ---")
EIRP_dBm = P_tx_dBm - L_filt_tx_dB - L_cable_tx_dB + G_tx_dB - L_ant_tx_dB
print(f"Effective Isotropic Radiated Power (EIRP) Calculation:")
print(f"EIRP = Tx Power - Tx Filter Loss - Tx Cable Loss + Tx Antenna Gain - Tx Antenna Loss")
print(f"EIRP = {P_tx_dBm:.1f} dBm - {L_filt_tx_dB:.1f} dB - {L_cable_tx_dB:.1f} dB + {G_tx_dB:.1f} dB - {L_ant_tx_dB:.1f} dB = {EIRP_dBm:.1f} dBm\n")

# --- Step 2: Calculate Free Space Path Loss (FSPL) ---
FSPL_dB = 20 * math.log10(d_m) + 20 * math.log10(freq_Hz) + 20 * math.log10(4 * math.pi / c)
print(f"Free Space Path Loss (FSPL) Calculation:")
print(f"FSPL = 20*log10(distance_m) + 20*log10(frequency_Hz) + 20*log10(4*pi/c)")
print(f"FSPL = 20*log10({d_m:.0f}) + 20*log10({freq_Hz:.0e}) - 147.55 = {FSPL_dB:.2f} dB\n")

# --- Step 3: Calculate Received Signal Power at LNA Input ---
# This is a convenient reference point for the SNR calculation.
P_received_at_ant_output_dBm = EIRP_dBm - FSPL_dB + G_rx_dB - L_ant_rx_dB
P_signal_at_LNA_input_dBm = P_received_at_ant_output_dBm - L_filt1_rx_dB
print(f"Received Signal Power Calculation (at LNA input):")
print(f"P_signal = EIRP - FSPL + Rx Ant Gain - Rx Ant Loss - Rx Filter Loss")
print(f"P_signal = {EIRP_dBm:.1f} dBm - {FSPL_dB:.2f} dB + {G_rx_dB:.1f} dB - {L_ant_rx_dB:.1f} dB - {L_filt1_rx_dB:.1f} dB = {P_signal_at_LNA_input_dBm:.2f} dBm\n")

# --- Step 4: Calculate Total Receiver Noise ---
print("--- 2. Noise Power Calculation ---")

# --- Step 4a: Calculate Thermal Noise Floor ---
Pn_watts = k * T * B_Hz
Pn_dBm = 10 * math.log10(Pn_watts / 0.001)
print(f"Thermal Noise Floor Calculation (at {T}K, {B_Hz/1000:.0f} kHz BW):")
print(f"Pn = 10 * log10(k * T * B / 0.001)")
print(f"Pn = 10 * log10({k:.2e} * {T} * {B_Hz:.0f} / 0.001) = {Pn_dBm:.2f} dBm\n")

# --- Step 4b: Calculate Receiver's Cascaded Noise Figure (from LNA onwards) ---
def from_db(x): return 10**(x / 10)
def to_db(x): return 10 * math.log10(x)

# Component parameters in linear scale needed for Friis formula
G_lna_lin = from_db(G_lna_dB)
F_lna_lin = from_db(NF_lna_dB)
G_mix_lin = from_db(-L_mix_dB) # Gain is reciprocal of loss
F_mix_lin = from_db(NF_mix_dB)
G_filt_if_lin = from_db(-L_filt_if_dB)
F_filt_if_lin = from_db(L_filt_if_dB)

# Friis formula for cascade starting at LNA. Later stages have negligible contribution.
F_cascade_lin = F_lna_lin + (F_mix_lin - 1) / G_lna_lin + (F_filt_if_lin - 1) / (G_lna_lin * G_mix_lin)
NF_cascade_dB = to_db(F_cascade_lin)

print(f"Receiver Noise Figure Calculation (from LNA onwards):")
print(f"NF_cascade = NF_lna + (NF_mixer - 1)/G_lna + ...")
print(f"Using values: NF_lna={NF_lna_dB:.1f} dB, G_lna={G_lna_dB:.1f} dB, NF_mixer={NF_mix_dB:.1f} dB")
print(f"Resulting cascaded Noise Figure from LNA onwards = {NF_cascade_dB:.2f} dB\n")

# --- Step 4c: Calculate Total Referred Noise Power at LNA Input ---
P_noise_total_at_LNA_input_dBm = Pn_dBm + NF_cascade_dB
print(f"Total Noise Power Calculation (at LNA input):")
print(f"P_noise = Thermal Noise Floor + Cascaded Noise Figure")
print(f"P_noise = {Pn_dBm:.2f} dBm + {NF_cascade_dB:.2f} dB = {P_noise_total_at_LNA_input_dBm:.2f} dBm\n")

# --- Step 5: Final SNR Calculation ---
print("--- 3. Final SNR Calculation ---")
SNR_dB = P_signal_at_LNA_input_dBm - P_noise_total_at_LNA_input_dBm
print(f"SNR = Signal Power at LNA input - Noise Power at LNA input")
print(f"SNR = {P_signal_at_LNA_input_dBm:.2f} dBm - ({P_noise_total_at_LNA_input_dBm:.2f} dBm)")
print(f"Resulting SNR = {SNR_dB:.2f} dB")
<<<28.27>>>