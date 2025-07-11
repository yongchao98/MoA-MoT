import math

# A function to convert dB values to linear scale for noise calculations
def db_to_linear(db):
    """Converts a value from decibels (dB) to a linear scale."""
    return 10**(db / 10)

# --- Given Parameters ---
# Transmitter (Tx)
p_tx_dbm = 30.0
g_tx_db = 20.0
l_tx_ant_db = 1.0
l_tx_filt_db = 1.0
l_tx_cable_db = 1.0

# Path and Signal
freq_ghz = 24.0
dist_km = 10.0
bw_khz = 100.0

# Receiver (Rx)
g_rx_db = 1.0
l_rx_ant_db = 0.5
l_rx_filt1_db = 1.0
g_lna_db = 36.0
nf_lna_db = 2.0
l_mix_db = 9.0  # Mixer loss
nf_mix_db = 9.0 # Mixer Noise Figure is assumed to be its conversion loss

# Environment
temp_k = 300.0
k_boltzmann = 1.380649e-23  # Boltzmann's constant in J/K

print("--- Link Budget and SNR Calculation ---")

# Step 1: Calculate Effective Isotropic Radiated Power (EIRP)
print("\nStep 1: Calculate EIRP (Effective Isotropic Radiated Power)")
l_tx_total_db = l_tx_ant_db + l_tx_filt_db + l_tx_cable_db
eirp_dbm = p_tx_dbm + g_tx_db - l_tx_total_db
print(f"Total Tx Loss = {l_tx_ant_db:.1f} dB (ant) + {l_tx_filt_db:.1f} dB (filter) + {l_tx_cable_db:.1f} dB (cable) = {l_tx_total_db:.2f} dB")
print(f"EIRP = Tx Power + Tx Gain - Total Tx Loss = {p_tx_dbm:.2f} dBm + {g_tx_db:.2f} dB - {l_tx_total_db:.2f} dB = {eirp_dbm:.2f} dBm")

# Step 2: Calculate Free Space Path Loss (FSPL)
print("\nStep 2: Calculate Free Space Path Loss (FSPL)")
freq_mhz = freq_ghz * 1000
fspl_db = 32.44 + 20 * math.log10(freq_mhz) + 20 * math.log10(dist_km)
print(f"FSPL = 32.44 + 20*log10({freq_mhz:.0f} MHz) + 20*log10({dist_km:.1f} km) = {fspl_db:.2f} dB")

# Step 3: Calculate Received Signal Power (at Rx antenna terminals)
print("\nStep 3: Calculate Received Signal Power (S)")
rx_signal_power_dbm = eirp_dbm - fspl_db + g_rx_db
print(f"S = EIRP - FSPL + Rx Gain = {eirp_dbm:.2f} dBm - {fspl_db:.2f} dB + {g_rx_db:.2f} dB = {rx_signal_power_dbm:.2f} dBm")

# Step 4: Calculate System Noise Figure (NF_sys)
print("\nStep 4: Calculate System Noise Figure (NF_sys)")
g_lna_lin = db_to_linear(g_lna_db)
nf_lna_lin = db_to_linear(nf_lna_db)
nf_mix_lin = db_to_linear(nf_mix_db)
# Noise figure of the cascade from the LNA onwards (Friis formula)
nf_cascaded_lin = nf_lna_lin + (nf_mix_lin - 1) / g_lna_lin
nf_cascaded_db = 10 * math.log10(nf_cascaded_lin)
nf_sys_db = l_rx_ant_db + l_rx_filt1_db + nf_cascaded_db
print(f"NF_cascaded (from LNA) = 10*log10({nf_lna_lin:.2f} + ({nf_mix_lin:.2f} - 1) / {g_lna_lin:.0f}) = {nf_cascaded_db:.2f} dB")
print(f"NF_sys = L_rx_ant + L_rx_filt1 + NF_cascaded = {l_rx_ant_db:.2f} dB + {l_rx_filt1_db:.2f} dB + {nf_cascaded_db:.2f} dB = {nf_sys_db:.2f} dB")

# Step 5: Calculate Total Noise Power (N)
print("\nStep 5: Calculate Total Noise Power (N)")
bw_hz = bw_khz * 1000
thermal_noise_density_dbm_hz = 10 * math.log10(k_boltzmann * temp_k * 1000)
thermal_noise_dbm = thermal_noise_density_dbm_hz + 10 * math.log10(bw_hz)
total_noise_power_dbm = thermal_noise_dbm + nf_sys_db
print(f"Thermal Noise = 10*log10(k*T*BW*1000) = {thermal_noise_dbm:.2f} dBm")
print(f"N = Thermal Noise + NF_sys = {thermal_noise_dbm:.2f} dBm + {nf_sys_db:.2f} dB = {total_noise_power_dbm:.2f} dBm")

# Step 6: Calculate Final Signal-to-Noise Ratio (SNR)
print("\nStep 6: Calculate Final Signal-to-Noise Ratio (SNR)")
snr_db = rx_signal_power_dbm - total_noise_power_dbm
print(f"SNR = S - N")
print(f"SNR = {rx_signal_power_dbm:.2f} dBm - ({total_noise_power_dbm:.2f} dBm) = {snr_db:.2f} dB")

final_answer = snr_db
print(f"\nFinal calculated SNR is {final_answer:.2f} dB.")
print(f"<<<{final_answer:.2f}>>>")