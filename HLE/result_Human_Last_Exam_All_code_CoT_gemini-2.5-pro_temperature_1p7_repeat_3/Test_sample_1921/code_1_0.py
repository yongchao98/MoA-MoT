import math

# Step 0: Define constants and given values
P_tx_dBm = 30.0  # Transmitter Power in dBm
G_tx_dB = 20.0  # Tx Antenna Gain in dB
L_tx_ant_dB = 1.0  # Tx Antenna Loss in dB
L_tx_filter_dB = 1.0  # Tx Filter Loss in dB
L_tx_cable_dB = 1.0  # Tx Cable Loss in dB

freq_Hz = 24e9  # Frequency in Hz (24 GHz)
dist_m = 10e3  # Distance in meters (10 km)
c_light_mps = 299792458  # Speed of light in m/s

G_rx_dB = 1.0  # Rx Antenna Gain in dB
L_rx_ant_dB = 0.5  # Rx Antenna Loss in dB
L_rx_filter_in_dB = 1.0  # Rx Input Filter Loss in dB

LNA_G_dB = 36.0  # LNA Gain in dB
LNA_NF_dB = 2.0  # LNA Noise Figure in dB

Mixer_L_dB = 9.0  # Mixer Conversion Loss in dB
Mixer_NF_dB = 9.0 # Noise figure of a passive mixer is its loss

IF_filter_L_dB = 1.0 # IF Filter Loss
IF_filter_NF_dB = 1.0

IF_amp_G_dB = 23.0 # IF Amplifier Gain
IF_amp_NF_dB = 0.0 # Negligible NF

Out_filter_L_dB = 1.0 # Output Filter Loss
Out_filter_NF_dB = 1.0

T_K = 300.0  # Ambient Temperature in Kelvin
k_Boltzmann = 1.380649e-23  # Boltzmann's constant
BW_Hz = 100e3  # Signal Bandwidth in Hz (100 kHz)

# Helper functions for dB conversion
def db_to_linear(db_val):
    """Converts a dB value to a linear ratio."""
    return 10**(db_val / 10.0)

def linear_to_db(linear_val):
    """Converts a linear ratio to a dB value."""
    return 10 * math.log10(linear_val)

# Step 1: Calculate Transmitter EIRP (Effective Isotropic Radiated Power)
L_tx_total_dB = L_tx_ant_dB + L_tx_filter_dB + L_tx_cable_dB
EIRP_dBm = P_tx_dBm + G_tx_dB - L_tx_total_dB
print("--- 1. Signal Power Calculation ---")
print(f"Transmitter EIRP = P_tx + G_tx_ant - L_tx_total")
print(f"EIRP = {P_tx_dBm} dBm + {G_tx_dB} dB - ({L_tx_ant_dB} + {L_tx_filter_dB} + {L_tx_cable_dB}) dB = {EIRP_dBm:.2f} dBm\n")

# Step 2: Calculate Free Space Path Loss (FSPL)
fspl_dB = 20 * math.log10(dist_m) + 20 * math.log10(freq_Hz) + 20 * math.log10(4 * math.pi / c_light_mps)
print("Free Space Path Loss (FSPL) = 20*log10(distance_m) + 20*log10(frequency_Hz) + 20*log10(4*pi/c)")
print(f"FSPL = 20*log10({dist_m:.0f}) + 20*log10({freq_Hz:.0f}) + 20*log10(4*pi/{c_light_mps}) = {fspl_dB:.2f} dB\n")

# Step 3: Calculate Received Power at antenna port (Prx)
Prx_dBm = EIRP_dBm - fspl_dB + G_rx_dB
print("Received Power (Prx) = EIRP - FSPL + G_rx_ant")
print(f"Prx = {EIRP_dBm:.2f} dBm - {fspl_dB:.2f} dB + {G_rx_dB} dB = {Prx_dBm:.2f} dBm\n")

# Step 4: Calculate Total Receiver Noise Figure (NF_rx) using Friis formula
# The receiver chain is considered from the antenna port onwards.
# Stages: (Gain_dB, Noise_Figure_dB)
stages = [
    (-L_rx_ant_dB, L_rx_ant_dB),
    (-L_rx_filter_in_dB, L_rx_filter_in_dB),
    (LNA_G_dB, LNA_NF_dB),
    (-Mixer_L_dB, Mixer_NF_dB),
    (-IF_filter_L_dB, IF_filter_NF_dB),
    (IF_amp_G_dB, IF_amp_NF_dB),
    (-Out_filter_L_dB, Out_filter_NF_dB)
]

# Friis Formula: F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
G_cumulative_lin = 1.0
F_total_lin = 0.0

for i, (g_db, nf_db) in enumerate(stages):
    g_lin = db_to_linear(g_db)
    f_lin = db_to_linear(nf_db)
    
    if i == 0:
        F_total_lin += f_lin
    else:
        F_total_lin += (f_lin - 1) / G_cumulative_lin
    
    G_cumulative_lin *= g_lin

NF_total_dB = linear_to_db(F_total_lin)
print("--- 2. Noise Power Calculation ---")
print("Receiver Noise Figure (NF_rx) is calculated using the Friis formula for cascaded stages.")
print(f"The combined Noise Figure of the receiver components is: {NF_total_dB:.2f} dB\n")

# Step 5: Calculate Total Noise Power (N_total)
# Convert thermal noise power (k*T*B) to dBm
P_noise_thermal_W = k_Boltzmann * T_K * BW_Hz
P_noise_thermal_dBm = 10 * math.log10(P_noise_thermal_W / 0.001)

N_total_dBm = P_noise_thermal_dBm + NF_total_dB
print("Total Noise Power (N_total) = Thermal Noise + Receiver Noise Figure")
print(f"Thermal Noise in 100kHz BW (kTB) = 10*log10(1.38e-23 * {T_K:.0f} * {BW_Hz:.0f} * 1000) = {P_noise_thermal_dBm:.2f} dBm")
print(f"N_total = {P_noise_thermal_dBm:.2f} dBm + {NF_total_dB:.2f} dB = {N_total_dBm:.2f} dBm\n")

# Step 6: Calculate final SNR
SNR_dB = Prx_dBm - N_total_dBm
print("--- 3. Final SNR Calculation ---")
print("SNR = Received Power (Prx) - Total Noise Power (N_total)")
print(f"SNR = {Prx_dBm:.2f} dBm - ({N_total_dBm:.2f} dBm)")
print(f"Final SNR = {SNR_dB:.2f} dB")
print("<<<28.26>>>")