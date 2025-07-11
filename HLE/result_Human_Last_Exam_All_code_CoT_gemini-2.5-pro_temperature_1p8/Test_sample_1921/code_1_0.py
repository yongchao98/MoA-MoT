import math

# --- Helper functions ---
def db_to_linear(db):
    """Converts a dB value to a linear scale."""
    return 10**(db / 10)

def linear_to_db(linear):
    """Converts a linear value to a dB scale."""
    return 10 * math.log10(linear)

# --- Given Parameters ---
# Transmitter (Tx)
P_tx_dbm = 30.0     # Tx Power in dBm
G_tx_ant_db = 20.0  # Tx Antenna Gain in dB
L_tx_ant_db = 1.0   # Tx Antenna Loss in dB
L_tx_filter_db = 1.0# Tx Filter Loss in dB
L_tx_cable_db = 1.0 # Tx Cable Loss in dB

# Path
f_hz = 24e9         # Frequency in Hz (24 GHz)
d_m = 10e3          # Distance in meters (10 km)

# Receiver (Rx)
G_rx_ant_db = 1.0   # Rx Antenna Gain in dB
L_rx_ant_db = 0.5   # Rx Antenna Loss in dB
L_rx_filter1_db = 1.0# Rx Input Filter Loss in dB

# Receiver Cascade Components for Noise Figure Calculation
# Stage 1: LNA
G_lna_db = 36.0     # LNA Gain in dB
NF_lna_db = 2.0     # LNA Noise Figure in dB
# Stage 2: Mixer
L_mixer_db = 9.0    # Mixer Conversion Loss in dB (acts as negative gain)
NF_mixer_db = 9.0   # NF of passive mixer is its loss
# Stage 3: IF Filter 1
L_if_filt1_db = 1.0 # IF Filter Loss in dB
NF_if_filt1_db = 1.0# NF of passive filter is its loss
# Stage 4: IF Amplifier
G_if_amp_db = 23.0  # IF Amp Gain in dB
NF_if_amp_db = 0.0  # Negligible NF
# Stage 5: Output Filter
L_if_filt2_db = 1.0 # Output Filter Loss in dB
NF_if_filt2_db = 1.0# NF of passive filter is its loss

# System
B_hz = 100e3        # Bandwidth in Hz (100 kHz)
T_kelvin = 300.0    # Ambient Temperature in Kelvin

# Constants
k = 1.380649e-23    # Boltzmann's constant in J/K
c = 299792458.0     # Speed of light in m/s

# --- Calculations ---

# 1. Calculate EIRP (Effective Isotropic Radiated Power)
L_tx_total_db = L_tx_ant_db + L_tx_filter_db + L_tx_cable_db
EIRP_dbm = P_tx_dbm + G_tx_ant_db - L_tx_total_db

# 2. Calculate Free Space Path Loss (FSPL)
# FSPL (dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
fspl_db = 20 * math.log10(d_m) + 20 * math.log10(f_hz) + 20 * math.log10(4 * math.pi / c)

# 3. Calculate Received Signal Power (S) at the LNA input
L_rx_pre_lna_db = L_rx_ant_db + L_rx_filter1_db
S_at_lna_input_dbm = EIRP_dbm - fspl_db + G_rx_ant_db - L_rx_pre_lna_db

# 4. Calculate Thermal Noise Power (N_thermal) in the signal bandwidth
N_thermal_watts = k * T_kelvin * B_hz
N_thermal_dbm = 10 * math.log10(N_thermal_watts / 0.001)

# 5. Calculate Receiver Noise Figure (NF_rx) for the chain starting at the LNA
# Use Friis formula: F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
components = [
    {'G_db': G_lna_db,      'NF_db': NF_lna_db},
    {'G_db': -L_mixer_db,   'NF_db': NF_mixer_db},
    {'G_db': -L_if_filt1_db,'NF_db': NF_if_filt1_db},
    {'G_db': G_if_amp_db,   'NF_db': NF_if_amp_db},
    {'G_db': -L_if_filt2_db,'NF_db': NF_if_filt2_db},
]

F_total_linear = 0.0
G_cascade_linear = 1.0

for comp in components:
    F_comp_linear = db_to_linear(comp['NF_db'])
    G_comp_linear = db_to_linear(comp['G_db'])
    
    noise_contribution = (F_comp_linear - 1) / G_cascade_linear
    F_total_linear += noise_contribution
    
    G_cascade_linear *= G_comp_linear

# The first component's F is added differently, adjust F_total by adding 1.
F_lna_chain_linear = F_total_linear + 1 
NF_lna_chain_db = linear_to_db(F_lna_chain_linear)

# 6. Calculate Final SNR
# SNR = S - N = S - (N_thermal + NF)
final_snr_db = S_at_lna_input_dbm - (N_thermal_dbm + NF_lna_chain_db)

# --- Print Results ---
print("--- Link Budget Calculation ---")
print(f"EIRP = {EIRP_dbm:.2f} dBm")
print(f"Free Space Path Loss = {fspl_db:.2f} dB")
print("\n--- SNR Calculation at LNA Input ---")
print(f"The final SNR is calculated using the formula: SNR = S - (N_thermal + NF_rx)")
print("\nWhere:")
print(f"S (Signal Power)           = {S_at_lna_input_dbm:.2f} dBm")
print(f"N_thermal (Thermal Noise)    = {N_thermal_dbm:.2f} dBm")
print(f"NF_rx (Receiver Noise Figure)= {NF_lna_chain_db:.2f} dB")
print("\nFinal Equation:")
print(f"SNR = {S_at_lna_input_dbm:.2f} dBm - ({N_thermal_dbm:.2f} dBm + {NF_lna_chain_db:.2f} dB)")
print(f"SNR = {S_at_lna_input_dbm:.2f} dBm - ({(N_thermal_dbm + NF_lna_chain_db):.2f} dBm)")
print("\n------------------------------------")
print(f"Resulting SNR = {final_snr_db:.2f} dB")
print("------------------------------------")

print(f"\n<<<SNR = {S_at_lna_input_dbm:.2f} - ({N_thermal_dbm:.2f} + {NF_lna_chain_db:.2f}) = {final_snr_db:.2f} dB>>>")
<<<28.27>>>