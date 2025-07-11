import math

def db_to_linear(db):
    """Converts a value from dB to a linear ratio."""
    return 10**(db / 10)

def linear_to_db(lin):
    """Converts a linear ratio to a value in dB."""
    return 10 * math.log10(lin)

# 0. System Constants and Parameters
# ----------------------------------------------------------------
# Constants
k = 1.38e-23  # Boltzmann's constant in J/K
c = 299792458   # Speed of light in m/s
T_K = 300       # Ambient temperature in Kelvin

# Transmitter (Tx) parameters
P_tx_dBm = 30.0       # Tx Power (dBm)
G_ant_tx_dB = 20.0    # Tx Antenna Gain (dB)
L_ant_tx_dB = 1.0     # Tx Antenna Loss (dB)
L_filt_tx_dB = 1.0    # Tx Output Filter Loss (dB)
L_cable_tx_dB = 1.0   # Tx Cable Loss (dB)

# Path parameters
freq_Hz = 24e9        # Frequency (Hz)
dist_m = 10e3         # Distance (m)

# Receiver (Rx) parameters
G_ant_rx_dB = 1.0     # Rx Antenna Gain (dB)
BW_Hz = 100e3         # Signal Bandwidth (Hz)

# Receiver chain components for noise figure calculation
# For passive components, Loss (dB) = Noise Figure (dB)
# Loss is represented as a negative gain
rx_chain = [
    {'name': 'Rx Antenna Loss', 'gain_dB': -0.5, 'nf_dB': 0.5},
    {'name': 'Rx Input Filter', 'gain_dB': -1.0, 'nf_dB': 1.0},
    {'name': 'LNA',             'gain_dB': 36.0, 'nf_dB': 2.0},
    {'name': 'Mixer',           'gain_dB': -9.0, 'nf_dB': 9.0}, # For a passive mixer, NF_dB = Loss_dB
    {'name': 'IF Filter',       'gain_dB': -1.0, 'nf_dB': 1.0},
    {'name': 'IF Amplifier',    'gain_dB': 23.0, 'nf_dB': 0.0}, # Negligible NF
    {'name': 'Output Filter',   'gain_dB': -1.0, 'nf_dB': 1.0}
]

# 1. Calculate Effective Isotropic Radiated Power (EIRP)
# ----------------------------------------------------------------
EIRP_dBm = P_tx_dBm + G_ant_tx_dB - L_ant_tx_dB - L_filt_tx_dB - L_cable_tx_dB
print(f"Step 1: EIRP = {P_tx_dBm} dBm + {G_ant_tx_dB} dB - {L_ant_tx_dB} dB - {L_filt_tx_dB} dB - {L_cable_tx_dB} dB = {EIRP_dBm:.2f} dBm")

# 2. Calculate Free Space Path Loss (FSPL)
# ----------------------------------------------------------------
# Using formula: FSPL(dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
FSPL_dB = 20 * math.log10(dist_m) + 20 * math.log10(freq_Hz) + 20 * math.log10(4 * math.pi / c)
print(f"Step 2: Free Space Path Loss (FSPL) = {FSPL_dB:.2f} dB")

# 3. Calculate Received Signal Power (Pr) at LNA input
# ----------------------------------------------------------------
# Pr = EIRP - FSPL + G_rx - L_rx_pre_lna
L_pre_lna_dB = rx_chain[0]['nf_dB'] + rx_chain[1]['nf_dB']
Pr_dBm = EIRP_dBm - FSPL_dB + G_ant_rx_dB - L_pre_lna_dB
print(f"Step 3: Received Signal Power (Pr) = {EIRP_dBm:.2f} dBm - {FSPL_dB:.2f} dB + {G_ant_rx_dB} dB - {L_pre_lna_dB} dB = {Pr_dBm:.2f} dBm")

# 4. Calculate System Noise Figure (NF_sys) using Friis Formula
# ----------------------------------------------------------------
# F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
F_sys_lin = 0
G_cascade_lin = 1.0

for component in rx_chain:
    F_lin = db_to_linear(component['nf_dB'])
    G_lin = db_to_linear(component['gain_dB'])
    
    F_sys_lin += (F_lin - 1) / G_cascade_lin
    G_cascade_lin *= G_lin

# The formula above calculates the noise from each stage referred to the input.
# To get the total noise factor, we need to add F1 at the beginning.
F1 = db_to_linear(rx_chain[0]['nf_dB'])
F_sys_lin = db_to_linear(rx_chain[0]['nf_dB']) # Fsys = F1
G1_lin = db_to_linear(rx_chain[0]['gain_dB'])
G_cascade_lin = G1_lin
# Recalculate properly
# F_sys = F_1 + (F_2-1)/G_1 + (F_3-1)/(G_1*G_2) + ...
F_sys_lin_terms = []
G_cumulative_lin = 1.0
for i, component in enumerate(rx_chain):
    f_lin = db_to_linear(component['nf_dB'])
    g_lin = db_to_linear(component['gain_dB'])
    
    if i == 0:
        term = f_lin
    else:
        term = (f_lin - 1) / G_cumulative_lin
        
    F_sys_lin_terms.append(term)
    G_cumulative_lin *= g_lin

F_sys_lin = sum(F_sys_lin_terms)
NF_sys_dB = linear_to_db(F_sys_lin)
print(f"Step 4: System Noise Figure (NF_sys) = {NF_sys_dB:.2f} dB")

# 5. Calculate Total Noise Power (N)
# ----------------------------------------------------------------
# N = k * T * B * F_sys
# In dBm: N_dBm = 10*log10(k*T*1000) + 10*log10(B) + NF_sys_dB
N0_dBm_per_Hz = 10 * math.log10(k * T_K * 1000)
BW_dBHz = 10 * math.log10(BW_Hz)
N_dBm = N0_dBm_per_Hz + BW_dBHz + NF_sys_dB
print(f"Step 5: Noise Power (N) = {N0_dBm_per_Hz:.2f} dBm/Hz + {BW_dBHz:.2f} dB-Hz + {NF_sys_dB:.2f} dB = {N_dBm:.2f} dBm")

# 6. Calculate Final Signal-to-Noise Ratio (SNR)
# ----------------------------------------------------------------
SNR_dB = Pr_dBm - N_dBm
print("\n--- Final Calculation ---")
print(f"SNR (dB) = Received Signal Power (dBm) - Total Noise Power (dBm)")
print(f"SNR (dB) = {Pr_dBm:.2f} dBm - ({N_dBm:.2f} dBm)")
print(f"Final SNR = {SNR_dB:.2f} dB")

print(f"\n<<<{SNR_dB:.2f}>>>")