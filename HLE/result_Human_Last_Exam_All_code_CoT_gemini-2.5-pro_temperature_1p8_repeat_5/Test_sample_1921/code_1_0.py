import math

# Step 0: Define constants and system parameters
# --- Constants
k = 1.380649e-23  # Boltzmann's constant in J/K
T = 300            # Ambient temperature in Kelvin
c = 299792458      # Speed of light in m/s

# --- System Parameters
# Transmitter
p_tx_dbm = 30.0
g_tx_db = 20.0
l_ant_tx_db = 1.0
l_filt_tx_db = 1.0
l_cable_db = 1.0
freq_hz = 24e9
bw_hz = 100e3

# Path
dist_m = 10e3

# Receiver
g_rx_db = 1.0
l_ant_rx_db = 0.5
l_filt1_rx_db = 1.0
g_lna_db = 36.0
nf_lna_db = 2.0
l_mixer_db = 9.0
l_filt_if_db = 1.0
g_if_db = 23.0
l_filt_out_db = 1.0

# --- Calculations ---
print("--- SNR Calculation Breakdown ---\n")

# Step 1: Calculate Effective Isotropic Radiated Power (EIRP)
print("1. Effective Isotropic Radiated Power (EIRP) Calculation:")
g_tx_net_db = g_tx_db - l_ant_tx_db
tx_losses_db = l_filt_tx_db + l_cable_db
eirp_dbm = p_tx_dbm - tx_losses_db + g_tx_net_db
print(f"   EIRP (dBm) = P_tx (dBm) - Tx Losses (dB) + Net Tx Antenna Gain (dB)")
print(f"   EIRP (dBm) = {p_tx_dbm:.1f} - ({l_filt_tx_db:.1f} + {l_cable_db:.1f}) + ({g_tx_db:.1f} - {l_ant_tx_db:.1f}) = {eirp_dbm:.2f} dBm\n")

# Step 2: Calculate Free Space Path Loss (FSPL)
print("2. Free Space Path Loss (FSPL) Calculation:")
# Using the formula FSPL(dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
fspl_db = 20 * math.log10(dist_m) + 20 * math.log10(freq_hz) + 20 * math.log10(4 * math.pi / c)
print(f"   FSPL (dB) = 20*log10({dist_m:.0f}) + 20*log10({freq_hz:.0f}) + 20*log10(4*pi/c)")
print(f"   FSPL (dB) = {20 * math.log10(dist_m):.2f} + {20 * math.log10(freq_hz):.2f} - {abs(20 * math.log10(4 * math.pi / c)):.2f} = {fspl_db:.2f} dB\n")

# Step 3: Calculate Received Signal Power (S) at the receiver input
print("3. Received Signal Power (S) Calculation:")
g_rx_net_db = g_rx_db - l_ant_rx_db
s_dbm = eirp_dbm - fspl_db + g_rx_net_db
print(f"   Signal (dBm) = EIRP (dBm) - FSPL (dB) + Net Rx Antenna Gain (dB)")
print(f"   Signal (dBm) = {eirp_dbm:.2f} - {fspl_db:.2f} + ({g_rx_db:.1f} - {l_ant_rx_db:.1f}) = {s_dbm:.2f} dBm\n")

# Step 4: Calculate Total System Noise Figure (NF_sys) using Friis Formula
print("4. System Noise Figure (NF_sys) Calculation:")
# Reference point: Input to the 150 MHz receiver filter
# Chain: Filter1 -> LNA -> Mixer -> IF Filter -> IF Amp

# Convert dB to linear values
# A loss 'L' in dB corresponds to a gain of '-L' dB.
# The noise figure of a passive component equals its loss.
# Component 1: Input Filter
g1_lin = 10**(-l_filt1_rx_db / 10)
f1_lin = 10**(l_filt1_rx_db / 10)

# Component 2: LNA
g2_lin = 10**(g_lna_db / 10)
f2_lin = 10**(nf_lna_db / 10)

# Component 3: Mixer (passive mixer assumed, NF = Loss)
g3_lin = 10**(-l_mixer_db / 10)
f3_lin = 10**(l_mixer_db / 10)

# Subsequent components have negligible impact due to high gain of LNA, but we calculate them for completeness.
g4_lin = 10**(-l_filt_if_db / 10)
f4_lin = 10**(l_filt_if_db / 10)

g5_lin = 10**(g_if_db / 10)
f5_lin = 1 # Negligible NF = 0 dB

# Friis Formula: F_sys = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
gain_cascaded = 1
f_sys_lin = 0

term1 = f1_lin
f_sys_lin += term1
gain_cascaded *= g1_lin

term2 = (f2_lin - 1) / gain_cascaded
f_sys_lin += term2
gain_cascaded *= g2_lin

term3 = (f3_lin - 1) / gain_cascaded
f_sys_lin += term3
gain_cascaded *= g3_lin

term4 = (f4_lin - 1) / gain_cascaded
f_sys_lin += term4

nf_sys_db = 10 * math.log10(f_sys_lin)

print(f"   F_sys = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...")
print(f"   F_sys = {f1_lin:.3f} + ({f2_lin:.3f}-1)/{g1_lin:.3f} + ({f3_lin:.3f}-1)/({g1_lin:.3f}*{g2_lin:.1f}) + ...")
print(f"   F_sys = {term1:.3f} + {term2:.3f} + {term3:.3f} + ... = {f_sys_lin:.3f}")
print(f"   NF_sys (dB) = 10*log10({f_sys_lin:.3f}) = {nf_sys_db:.2f} dB\n")


# Step 5: Calculate Total Noise Power (N) in the signal bandwidth
print("5. Total Noise Power (N) Calculation:")
n_thermal_dbm_hz = 10 * math.log10(k * T) + 30 # convert to dBm
n_dbm = n_thermal_dbm_hz + 10 * math.log10(bw_hz) + nf_sys_db
print(f"   Noise (dBm) = Thermal Noise (dBm/Hz) + 10*log10(Bandwidth) + NF_sys (dB)")
print(f"   Noise (dBm) = {n_thermal_dbm_hz:.2f} + {10*math.log10(bw_hz):.1f} + {nf_sys_db:.2f} = {n_dbm:.2f} dBm\n")

# Step 6: Calculate final SNR
print("6. Final Signal-to-Noise Ratio (SNR) Calculation:")
snr_db = s_dbm - n_dbm
print(f"   SNR (dB) = Signal (dBm) - Noise (dBm)")
print(f"   SNR (dB) = {s_dbm:.2f} - ({n_dbm:.2f}) = {snr_db:.2f} dB\n")

# Final Answer
print(f"The resulting SNR of the received signal is {snr_db:.2f} dB.")
final_snr = round(snr_db, 2)