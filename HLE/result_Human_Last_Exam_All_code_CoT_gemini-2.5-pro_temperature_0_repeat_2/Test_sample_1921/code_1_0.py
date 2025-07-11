import math

def db_to_linear(db):
    """Converts a dB value to a linear ratio."""
    return 10**(db / 10)

def linear_to_db(linear):
    """Converts a linear ratio to a dB value."""
    return 10 * math.log10(linear)

# --- 0. System Parameters ---
# Transmitter
p_tx_dbm = 30.0      # Transmitter power in dBm
g_tx_db = 20.0       # Tx antenna gain in dB
l_tx_ant_db = 1.0    # Tx antenna loss in dB
l_tx_filt_db = 1.0   # Tx filter loss in dB
l_tx_cable_db = 1.0  # Tx cable loss in dB

# Path
freq_ghz = 24.0      # Frequency in GHz
dist_km = 10.0       # Distance in km
bw_khz = 100.0       # Signal bandwidth in kHz

# Receiver
g_rx_db = 1.0        # Rx antenna gain in dB
l_rx_ant_db = 0.5    # Rx antenna loss in dB
l_rx_filt1_db = 1.0  # Rx input filter loss in dB
g_lna_db = 36.0      # LNA gain in dB
nf_lna_db = 2.0      # LNA Noise Figure in dB
l_mixer_db = 9.0     # Mixer conversion loss in dB (negative gain)
nf_mixer_db = 9.0    # Mixer Noise Figure (assumed equal to loss)
l_if_filt_db = 1.0   # IF filter loss in dB
g_if_amp_db = 23.0   # IF amplifier gain in dB
nf_if_amp_db = 0.0   # IF amplifier Noise Figure (negligible)
l_out_filt_db = 1.0  # Output filter loss in dB

# Constants
T_kelvin = 300.0     # Ambient temperature in Kelvin
k_boltzmann = 1.380649e-23 # Boltzmann's constant in J/K

# --- 1. Calculate EIRP (Effective Isotropic Radiated Power) ---
l_tx_total_db = l_tx_ant_db + l_tx_filt_db + l_tx_cable_db
eirp_dbm = p_tx_dbm + g_tx_db - l_tx_total_db
print("Step 1: Calculate Transmitter EIRP")
print(f"  Tx Power ({p_tx_dbm:.1f} dBm) + Tx Gain ({g_tx_db:.1f} dB) - Total Tx Loss ({l_tx_total_db:.1f} dB) = {eirp_dbm:.2f} dBm\n")

# --- 2. Calculate Free Space Path Loss (FSPL) ---
# FSPL (dB) = 20*log10(d) + 20*log10(f) + 20*log10(4π/c)
# Simplified for km and GHz: FSPL (dB) = 20*log10(d_km) + 20*log10(f_ghz) + 92.45
fspl_db = 20 * math.log10(dist_km) + 20 * math.log10(freq_ghz) + 92.45
print("Step 2: Calculate Free Space Path Loss (FSPL)")
print(f"  FSPL for {freq_ghz:.1f} GHz at {dist_km:.1f} km = {fspl_db:.2f} dB\n")

# --- 3. Calculate Received Signal Power (at receiver input reference point) ---
# Reference point is after Rx antenna gain, before any Rx losses.
p_rx_signal_dbm = eirp_dbm - fspl_db + g_rx_db
print("Step 3: Calculate Received Signal Power (P_rx)")
print(f"  EIRP ({eirp_dbm:.2f} dBm) - FSPL ({fspl_db:.2f} dB) + Rx Gain ({g_rx_db:.1f} dB) = {p_rx_signal_dbm:.2f} dBm\n")

# --- 4. Calculate System Noise Figure (NF_sys) ---
# Using Friis formula for cascaded noise figure: F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
# Stages are defined from the reference point onwards.
stages = [
    {'name': 'Rx Ant+Filt Loss', 'gain_db': -(l_rx_ant_db + l_rx_filt1_db), 'nf_db': (l_rx_ant_db + l_rx_filt1_db)},
    {'name': 'LNA',              'gain_db': g_lna_db, 'nf_db': nf_lna_db},
    {'name': 'Mixer',            'gain_db': -l_mixer_db, 'nf_db': nf_mixer_db},
    {'name': 'IF Filter',        'gain_db': -l_if_filt_db, 'nf_db': l_if_filt_db},
    {'name': 'IF Amp',           'gain_db': g_if_amp_db, 'nf_db': nf_if_amp_db},
    {'name': 'Output Filter',    'gain_db': -l_out_filt_db, 'nf_db': l_out_filt_db}
]

total_noise_factor = 0
cascaded_gain_linear = 1.0

for i, stage in enumerate(stages):
    gain_linear = db_to_linear(stage['gain_db'])
    noise_factor_linear = db_to_linear(stage['nf_db'])
    
    if i == 0:
        stage_contribution = noise_factor_linear
    else:
        stage_contribution = (noise_factor_linear - 1) / cascaded_gain_linear
        
    total_noise_factor += stage_contribution
    cascaded_gain_linear *= gain_linear

nf_sys_db = linear_to_db(total_noise_factor)
print("Step 4: Calculate System Noise Figure (NF_sys)")
print(f"  The total noise figure of the cascaded receiver components is {nf_sys_db:.2f} dB\n")

# --- 5. Calculate Total Noise Power ---
bw_hz = bw_khz * 1000
# Thermal noise in dBm = 10*log10(k * T * B * 1000)
thermal_noise_dbm = 10 * math.log10(k_boltzmann * T_kelvin * bw_hz * 1000)
total_noise_dbm = thermal_noise_dbm + nf_sys_db
print("Step 5: Calculate Total Noise Power (N)")
print(f"  Thermal Noise in {bw_khz:.0f} kHz BW ({thermal_noise_dbm:.2f} dBm) + System NF ({nf_sys_db:.2f} dB) = {total_noise_dbm:.2f} dBm\n")

# --- 6. Calculate Final SNR ---
snr_db = p_rx_signal_dbm - total_noise_dbm
print("Step 6: Calculate Final Signal-to-Noise Ratio (SNR)")
print(f"  The final equation is: SNR = P_rx - N")
print(f"  SNR = {p_rx_signal_dbm:.2f} dBm - ({total_noise_dbm:.2f} dBm)")
print(f"  Resulting SNR = {snr_db:.2f} dB")

# --- Final Answer ---
# The final answer is the SNR value.
final_answer = f"<<<{snr_db:.2f}>>>"
print(final_answer)