import math

def db_to_linear(db):
    """Converts a dB value to a linear ratio."""
    return 10**(db / 10)

def linear_to_db(linear):
    """Converts a linear ratio to a dB value."""
    return 10 * math.log10(linear)

# --- GIVEN PARAMETERS ---
# Transmitter (Tx)
P_tx_dBm = 30.0
G_tx_ant_dB = 20.0
L_tx_ant_dB = 1.0
L_tx_filter_dB = 1.0
L_tx_cable_dB = 1.0
f_GHz = 24.0
BW_kHz = 100.0

# Path
d_km = 10.0

# Receiver (Rx)
G_rx_ant_dB = 1.0
L_rx_ant_dB = 0.5
L_rx_in_filter_dB = 1.0
G_lna_dB = 36.0
NF_lna_dB = 2.0
L_mixer_dB = 9.0  # Conversion Loss
NF_mixer_dB = 9.0 # Assume NF equals conversion loss
L_if_filter_dB = 1.0
G_if_amp_dB = 23.0
NF_if_amp_dB = 0.0 # Negligible
L_out_filter_dB = 1.0 # This component's NF contribution is negligible after high gain stages

# System
T_K = 300.0

# Constants
k = 1.380649e-23  # Boltzmann's constant in J/K
c = 299792458     # Speed of light in m/s

# --- STEP 1: Calculate Effective Isotropic Radiated Power (EIRP) ---
print("--- Step 1: Transmitter Power Calculation ---")
L_tx_total_dB = L_tx_ant_dB + L_tx_filter_dB + L_tx_cable_dB
EIRP_dBm = P_tx_dBm + G_tx_ant_dB - L_tx_total_dB
print(f"Tx Power: {P_tx_dBm} dBm")
print(f"Tx Antenna Gain: {G_tx_ant_dB} dB")
print(f"Total Tx Losses (antenna + filter + cable): {L_tx_total_dB:.1f} dB")
print(f"EIRP = {P_tx_dBm:.1f} dBm + {G_tx_ant_dB:.1f} dB - {L_tx_total_dB:.1f} dB = {EIRP_dBm:.2f} dBm\n")

# --- STEP 2: Calculate Free Space Path Loss (FSPL) ---
print("--- Step 2: Free Space Path Loss (FSPL) Calculation ---")
d_m = d_km * 1000
f_hz = f_GHz * 1e9
fspl_dB = 20 * math.log10(d_m) + 20 * math.log10(f_hz) + 20 * math.log10(4 * math.pi / c)
print(f"Distance: {d_km} km")
print(f"Frequency: {f_GHz} GHz")
print(f"FSPL = {fspl_dB:.2f} dB\n")

# --- STEP 3: Calculate Received Signal Power (S) ---
print("--- Step 3: Received Signal Power (S) Calculation ---")
# Signal power at the input of the Rx chain (before the first filter)
effective_G_rx_ant_dB = G_rx_ant_dB - L_rx_ant_dB
S_at_rx_input_dBm = EIRP_dBm - fspl_dB + effective_G_rx_ant_dB
print("Signal Power is calculated at the input to the receiver's first electronic component.")
print(f"S (dBm) = EIRP - FSPL + Rx_Antenna_Gain_Net")
print(f"S (dBm) = {EIRP_dBm:.2f} dBm - {fspl_dB:.2f} dB + {effective_G_rx_ant_dB:.1f} dB = {S_at_rx_input_dBm:.2f} dBm\n")

# --- STEP 4: Calculate Receiver Total Noise Figure (NF) ---
print("--- Step 4: Receiver Total Noise Figure (NF) Calculation ---")
# The cascade starts at the input filter. Noise is referred to this point.
# Component 1: Input Filter
G1_lin = db_to_linear(-L_rx_in_filter_dB)
F1_lin = db_to_linear(L_rx_in_filter_dB) # NF of a passive component equals its loss

# Component 2: LNA
G2_lin = db_to_linear(G_lna_dB)
F2_lin = db_to_linear(NF_lna_dB)

# Component 3: Mixer
G3_lin = db_to_linear(-L_mixer_dB)
F3_lin = db_to_linear(NF_mixer_dB)

# Friis Formula for Noise Factor: F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
term1 = F1_lin
term2 = (F2_lin - 1) / G1_lin
term3 = (F3_lin - 1) / (G1_lin * G2_lin)
# Subsequent terms are negligible due to high gain of LNA and will not be shown.

F_total_lin = term1 + term2 + term3
NF_total_dB = linear_to_db(F_total_lin)
print("Using Friis formula for cascaded noise figure: F_total = F1 + (F2-1)/G1 + ...")
print(f"Term 1 (Input Filter): {term1:.3f}")
print(f"Term 2 (LNA): {term2:.3f}")
print(f"Term 3 (Mixer): {term3:.3f}")
print(f"Total Noise Factor (linear) F_total = {F_total_lin:.3f}")
print(f"Total Noise Figure NF_total = {NF_total_dB:.2f} dB\n")

# --- STEP 5: Calculate Total Noise Power (N) ---
print("--- Step 5: Total Noise Power (N) Calculation ---")
BW_hz = BW_kHz * 1000
# Thermal Noise Power N_thermal = k * T * B
N_thermal_dBm = linear_to_db((k * T_K * BW_hz) * 1000)
# Total Noise Power N = N_thermal + NF_total
N_dBm = N_thermal_dBm + NF_total_dB
print(f"Noise is calculated in the signal bandwidth of {BW_kHz} kHz at {T_K}K.")
print(f"Thermal Noise Floor = {N_thermal_dBm:.2f} dBm")
print(f"Total Noise Power (N) = Thermal_Noise_Floor + NF_total")
print(f"N (dBm) = {N_thermal_dBm:.2f} dBm + {NF_total_dB:.2f} dB = {N_dBm:.2f} dBm\n")

# --- STEP 6: Calculate Signal-to-Noise Ratio (SNR) ---
print("--- Step 6: Final Signal-to-Noise Ratio (SNR) Calculation ---")
SNR_dB = S_at_rx_input_dBm - N_dBm
print(f"SNR (dB) = Signal Power (dBm) - Noise Power (dBm)")
print(f"SNR (dB) = {S_at_rx_input_dBm:.2f} - ({N_dBm:.2f})")
print(f"Final SNR = {SNR_dB:.2f} dB")
print(f"\n<<<{SNR_dB:.2f}>>>")