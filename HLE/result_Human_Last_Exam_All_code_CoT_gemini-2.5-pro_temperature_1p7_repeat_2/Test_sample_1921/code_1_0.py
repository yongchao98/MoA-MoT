import math

# This script calculates the final Signal-to-Noise Ratio (SNR) for the given microwave link.

# --- Constants and Given Values ---
# Constants
k = 1.380649e-23  # Boltzmann's constant in J/K
c = 299792458      # Speed of light in m/s
T = 300            # Ambient temperature in Kelvin

# Transmitter
P_tx_dbm = 30.0
G_tx_db = 20.0
L_tx_total_db = 1.0 + 1.0 + 1.0 # Antenna loss + filter loss + cable loss

# Path
freq_hz = 24e9
distance_m = 10e3
signal_bw_hz = 100e3

# Receiver
G_rx_db = 1.0
# Receiver chain components for noise figure calculation
rx_components = [
    {'name': 'Rx Ant Loss', 'loss_db': 0.5},
    {'name': 'Rx Input Filter', 'loss_db': 1.0},
    {'name': 'LNA', 'gain_db': 36.0, 'nf_db': 2.0},
    {'name': 'Mixer', 'loss_db': 9.0},
    {'name': 'IF Filter', 'loss_db': 1.0},
    {'name': 'IF Amp', 'gain_db': 23.0, 'nf_db': 0.0}, # Negligible NF
    {'name': 'Output Filter', 'loss_db': 1.0},
]

# --- Helper Functions ---
def db_to_linear(db_val):
    """Converts a dB value to a linear ratio."""
    return 10**(db_val / 10)

def linear_to_db(lin_val):
    """Converts a linear ratio to a dB value."""
    return 10 * math.log10(lin_val)

# --- Calculation Steps ---

# 1. Calculate Effective Isotropic Radiated Power (EIRP)
eirp_dbm = P_tx_dbm + G_tx_db - L_tx_total_db

# 2. Calculate Free Space Path Loss (FSPL)
# Using the formula: FSPL (dB) = 20*log10(d) + 20*log10(f) + 20*log10(4π/c)
fspl_db = 20 * math.log10(distance_m) + 20 * math.log10(freq_hz) + 20 * math.log10(4 * math.pi / c)

# 3. Calculate Received Signal Power (S_in) at the receiver's input reference point
S_in_dbm = eirp_dbm - fspl_db + G_rx_db

# 4. Calculate Thermal Noise Power (N_in) in the signal bandwidth
N_in_watts = k * T * signal_bw_hz
N_in_dbm = 10 * math.log10(N_in_watts / 0.001)

# 5. Calculate total receiver Noise Figure (NF_total) using Friis formula
F_total_lin = 0
G_cascade_lin = 1.0

for component in rx_components:
    if 'loss_db' in component:
        L_lin = db_to_linear(component['loss_db'])
        F_lin = L_lin # Noise Figure of a passive component is its loss
        G_lin = 1 / L_lin
    else: # It's an active component like an amplifier
        G_lin = db_to_linear(component['gain_db'])
        F_lin = db_to_linear(component['nf_db'])

    if F_total_lin == 0: # First component
        F_total_lin = F_lin
    else:
        F_total_lin += (F_lin - 1) / G_cascade_lin
    
    G_cascade_lin *= G_lin

NF_total_db = linear_to_db(F_total_lin)

# 6. Calculate the final SNR
# SNR (dB) = S_in_dbm - N_in_dbm - NF_total_db
snr_db = S_in_dbm - N_in_dbm - NF_total_db
snr_input_db = S_in_dbm - N_in_dbm

# --- Output the results ---
print("--- Link Budget Calculation ---")
print(f"Signal Power at Receiver Input (S_in) = {S_in_dbm:.2f} dBm")
print(f"Thermal Noise Power in Bandwidth (N_in) = {N_in_dbm:.2f} dBm")
print(f"Total Receiver Noise Figure (NF_total) = {NF_total_db:.2f} dB")
print("\n--- Final SNR Equation ---")
print("Final SNR = (S_in - N_in) - NF_total")
print(f"Final SNR = ({S_in_dbm:.2f} dBm - ({N_in_dbm:.2f} dBm)) - {NF_total_db:.2f} dB")
print(f"Final SNR = {snr_input_db:.2f} dB - {NF_total_db:.2f} dB")
print(f"Resulting SNR = {snr_db:.2f} dB")