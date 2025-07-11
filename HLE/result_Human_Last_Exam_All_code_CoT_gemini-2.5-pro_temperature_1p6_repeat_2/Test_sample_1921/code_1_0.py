import math

def dB_to_linear(db):
    """Converts a dB value to a linear ratio."""
    return 10**(db / 10)

def linear_to_dB(linear):
    """Converts a linear ratio to a dB value."""
    return 10 * math.log10(linear)

# --- 1. System Parameters ---
# Transmitter
P_tx_dBm = 30
G_tx_dB = 20
L_tx_ant_dB = 1
L_tx_filt_dB = 1
L_cable_dB = 1

# Path
f_Hz = 24e9  # 24 GHz
d_km = 10
d_m = d_km * 1000

# Receiver
G_rx_dB = 1
B_Hz = 100e3 # 100 kHz signal bandwidth

# Receiver Chain for Noise Figure Calculation (in order)
# Component: (Gain in dB, Noise Figure in dB)
# For passive components (losses), Gain is -Loss and NF is Loss.
rx_chain = [
    ('Rx Ant Loss',  {'gain_dB': -0.5, 'nf_dB': 0.5}),
    ('Rx Filt Loss', {'gain_dB': -1.0, 'nf_dB': 1.0}),
    ('LNA',          {'gain_dB': 36,   'nf_dB': 2.0}),
    ('Mixer',        {'gain_dB': -9.0, 'nf_dB': 9.0}),
    ('IF Filt Loss', {'gain_dB': -1.0, 'nf_dB': 1.0}),
    ('IF Amp',       {'gain_dB': 23,   'nf_dB': 0.0}), # Negligible NF
]

# Environment
T_Kelvin = 300
k_Boltzmann = 1.380649e-23

# --- 2. Signal Power Calculation ---
print("--- Calculating Signal Power ---")

# a) Calculate Effective Isotropic Radiated Power (EIRP)
L_tx_total_dB = L_tx_ant_dB + L_tx_filt_dB + L_cable_dB
eirp_dBm = P_tx_dBm + G_tx_dB - L_tx_total_dB
print(f"Transmitter Power: {P_tx_dBm} dBm")
print(f"Transmit Antenna Gain: {G_tx_dB} dB")
print(f"Total Transmitter Losses: {L_tx_total_dB} dB")
print(f"Effective Isotropic Radiated Power (EIRP): {eirp_dBm:.2f} dBm")
print("-" * 20)

# b) Calculate Free Space Path Loss (FSPL)
# FSPL (dB) = 20*log10(d_m) + 20*log10(f_Hz) - 147.55
fspl_dB = 20 * math.log10(d_m) + 20 * math.log10(f_Hz) - 147.55
print(f"Distance: {d_km} km")
print(f"Frequency: {f_Hz / 1e9} GHz")
print(f"Free Space Path Loss (FSPL): {fspl_dB:.2f} dB")
print("-" * 20)

# c) Calculate Received Signal Power (S) at the antenna terminals
S_dBm = eirp_dBm - fspl_dB + G_rx_dB
print(f"Receiver Antenna Gain: {G_rx_dB} dB")
print(f"Received Signal Power (S): {eirp_dBm:.2f} dBm - {fspl_dB:.2f} dB + {G_rx_dB} dB = {S_dBm:.2f} dBm")
print("\n" + "="*35 + "\n")

# --- 3. Noise Power Calculation ---
print("--- Calculating Noise Power ---")

# a) Calculate System Noise Figure (NF_sys) using Friis formula
# F_sys = F_1 + (F_2-1)/G_1 + (F_3-1)/(G_1*G_2) + ...
F_sys_linear = 0
G_accum_linear = 1

print("Calculating System Noise Figure (NF_sys) from receiver chain:")
for name, params in rx_chain:
    F_i = dB_to_linear(params['nf_dB'])
    G_i = dB_to_linear(params['gain_dB'])
    
    contribution = (F_i - 1) / G_accum_linear
    F_sys_linear += contribution
    G_accum_linear *= G_i
    # The first component's full noise factor is added directly
    if name == 'Rx Ant Loss':
        F_sys_linear = F_i

NF_sys_dB = linear_to_dB(F_sys_linear)
print(f"Total System Noise Factor (F_sys): {F_sys_linear:.2f}")
print(f"Total System Noise Figure (NF_sys): {NF_sys_dB:.2f} dB")
print("-" * 20)

# b) Calculate thermal noise floor in the signal bandwidth
N_floor_watts = k_Boltzmann * T_Kelvin * B_Hz
N_floor_dBm = 10 * math.log10(N_floor_watts / 0.001)
print(f"Signal Bandwidth: {B_Hz/1e3} kHz")
print(f"Thermal Noise Floor in Bandwidth (k*T*B): {N_floor_dBm:.2f} dBm")
print("-" * 20)

# c) Calculate total noise power (N) at the same point as the signal
N_dBm = N_floor_dBm + NF_sys_dB
print(f"Total Noise Power (N): {N_floor_dBm:.2f} dBm + {NF_sys_dB:.2f} dB = {N_dBm:.2f} dBm")
print("\n" + "="*35 + "\n")

# --- 4. Final SNR Calculation ---
print("--- Calculating Final SNR ---")
SNR_dB = S_dBm - N_dBm
print(f"Final SNR (dB) = Signal Power (dBm) - Noise Power (dBm)")
print(f"Final SNR (dB) = {S_dBm:.2f} dBm - ({N_dBm:.2f} dBm) = {SNR_dB:.2f} dB")

# --- Final Answer ---
print(f"\n<<<The resulting SNR is {SNR_dB:.2f} dB>>>")
