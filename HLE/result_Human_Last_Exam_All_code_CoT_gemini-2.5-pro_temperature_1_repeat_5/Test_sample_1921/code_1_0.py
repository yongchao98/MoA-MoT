import math

# This script calculates the Signal-to-Noise Ratio (SNR) for the given microwave link.

# Plan:
# 1.  Calculate the transmitter's Effective Isotropic Radiated Power (EIRP).
#     EIRP = Tx Power + Antenna Gain - All Tx-side losses.
# 2.  Calculate the Free Space Path Loss (FSPL) for the given distance and frequency.
#     FSPL (dB) = 20*log10(distance) + 20*log10(frequency) + 20*log10(4*pi/c).
# 3.  Calculate the signal power arriving at the receiver's input terminals.
#     Received Power = EIRP - FSPL + Rx Antenna Gain - Rx Antenna Loss.
# 4.  Calculate the receiver's system noise figure (NF_sys) using the Friis formula.
# 5.  Calculate the thermal noise power (N) based on the signal bandwidth and ambient temperature.
#     N = k * T * B, where k is Boltzmann's constant.
# 6.  Calculate the total input-referred noise power by adding the system noise figure to the thermal noise power.
# 7.  Calculate the SNR by subtracting the total noise power from the received signal power.

# --- Constants and Parameters ---
# Constants
T = 300.0  # K (Ambient Temperature)
k = 1.38e-23  # J/K (Boltzmann's constant)
c = 299792458.0  # m/s (Speed of light)

# Transmitter Parameters
P_tx_dBm = 30.0  # dBm
G_tx_dB = 20.0  # dB
L_ant_tx_dB = 1.0  # dB
L_filter_tx_dB = 1.0  # dB
L_cable_tx_dB = 1.0  # dB
f_hz = 24e9  # Hz
B_hz = 100e3  # Hz

# Path Parameters
d_m = 10e3  # meters

# Receiver Parameters
G_rx_dB = 1.0  # dB
L_ant_rx_dB = 0.5  # dB
L_filter1_rx_dB = 1.0  # dB
G_lna_dB = 36.0  # dB
NF_lna_dB = 2.0  # dB
L_mixer_dB = 9.0  # dB
NF_mixer_dB = 9.0  # dB (Assuming NF = Conversion Loss for passive mixer)
# Subsequent components have negligible noise contribution due to high preceding gain

# --- Step 1: Calculate Signal Power ---
print("### Calculating Signal Power ###\n")

# EIRP Calculation
total_tx_loss_dB = L_ant_tx_dB + L_filter_tx_dB + L_cable_tx_dB
eirp_dbm = P_tx_dBm + G_tx_dB - total_tx_loss_dB
print("1. Effective Isotropic Radiated Power (EIRP):")
print(f"   EIRP = {P_tx_dBm} dBm (Tx Power) + {G_tx_dB} dB (Tx Ant Gain) - {L_filter_tx_dB} dB (Tx Filter Loss) - {L_cable_tx_dB} dB (Tx Cable Loss) - {L_ant_tx_dB} dB (Tx Ant Loss)")
print(f"   EIRP = {eirp_dbm:.2f} dBm\n")

# FSPL Calculation
fspl_db = 20 * math.log10(d_m) + 20 * math.log10(f_hz) + 20 * math.log10((4 * math.pi) / c)
print(f"2. Free Space Path Loss (FSPL) at {f_hz/1e9} GHz over {d_m/1e3} km:")
print(f"   FSPL = 20*log10({d_m:.0f} m) + 20*log10({f_hz:.0f} Hz) + 20*log10(4*pi/c)")
print(f"   FSPL = {fspl_db:.2f} dB\n")

# Received Signal Power Calculation
effective_rx_gain_db = G_rx_dB - L_ant_rx_dB
S_at_rx_input_dBm = eirp_dbm - fspl_db + effective_rx_gain_db
print(f"3. Received Signal Power (S) at Receiver Input:")
print(f"   S = {eirp_dbm:.2f} dBm (EIRP) - {fspl_db:.2f} dB (FSPL) + {G_rx_dB} dB (Rx Ant Gain) - {L_ant_rx_dB} dB (Rx Ant Loss)")
print(f"   S = {S_at_rx_input_dBm:.2f} dBm\n")


# --- Step 2: Calculate Noise Power ---
print("### Calculating Noise Power ###\n")

# Thermal Noise Calculation
thermal_noise_W = k * T * B_hz
thermal_noise_dBm = 10 * math.log10(thermal_noise_W / 0.001)
print(f"1. Thermal Noise Power (kTB) in {B_hz/1e3} kHz bandwidth at {T}K:")
print(f"   kTB = {k:.2e} J/K * {T} K * {B_hz:.0f} Hz")
print(f"   Thermal Noise = {thermal_noise_dBm:.2f} dBm\n")

# Receiver System Noise Figure Calculation
# Stage 1 (pre-LNA): Rx Antenna Loss + Rx Input Filter Loss
L_pre_lna_db = L_ant_rx_dB + L_filter1_rx_dB
F1_lin = 10**(L_pre_lna_db / 10) # NF of passive component = Loss
G1_lin = 10**(-L_pre_lna_db / 10)

# Stage 2: LNA
F2_lin = 10**(NF_lna_dB / 10)
G2_lin = 10**(G_lna_dB / 10)

# Stage 3: Mixer
F3_lin = 10**(NF_mixer_dB / 10)
G3_lin = 10**(-L_mixer_dB / 10)

# Friis Formula: F_sys = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
# Contributions from stages after the mixer are negligible due to the high gain of the LNA.
F_sys_lin = F1_lin + (F2_lin - 1) / G1_lin + (F3_lin - 1) / (G1_lin * G2_lin)
NF_sys_db = 10 * math.log10(F_sys_lin)

print(f"2. Receiver System Noise Figure (NF_sys):")
print(f"   Calculated using Friis formula for cascaded components.")
print(f"   - Loss before LNA = {L_ant_rx_dB} dB (Rx Ant Loss) + {L_filter1_rx_dB} dB (Rx Filter Loss) = {L_pre_lna_db:.1f} dB")
print(f"   - LNA Noise Figure = {NF_lna_dB} dB")
print(f"   - Mixer Noise Figure = {NF_mixer_dB} dB")
print(f"   Cascaded System Noise Figure (NF_sys) = {NF_sys_db:.2f} dB\n")

# Total Input-Referred Noise Power Calculation
N_total_dBm = thermal_noise_dBm + NF_sys_db
print(f"3. Total Input-Referred Noise Power (N):")
print(f"   N = {thermal_noise_dBm:.2f} dBm (Thermal Noise) + {NF_sys_db:.2f} dB (NF_sys)")
print(f"   N = {N_total_dBm:.2f} dBm\n")


# --- Step 3: Calculate Final SNR ---
print("### Calculating Final SNR ###\n")

SNR_dB = S_at_rx_input_dBm - N_total_dBm
print("Final Signal-to-Noise Ratio (SNR) Equation:")
print(f"SNR = Signal Power - Noise Power")
print(f"SNR = {S_at_rx_input_dBm:.2f} dBm - ({N_total_dBm:.2f} dBm)")
print(f"SNR = {SNR_dB:.2f} dB")
