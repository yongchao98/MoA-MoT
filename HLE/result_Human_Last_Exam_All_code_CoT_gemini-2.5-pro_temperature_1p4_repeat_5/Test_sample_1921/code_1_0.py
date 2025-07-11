import math

# --- 0. Constants and Helper Functions ---
# Boltzmann's constant in J/K
K_BOLTZMANN = 1.380649e-23
# Speed of light in m/s
C_LIGHT = 299792458

def db_to_linear(db):
    """Converts a dB value to a linear scale."""
    return 10**(db / 10)

def linear_to_db(linear):
    """Converts a linear value to a dB scale."""
    return 10 * math.log10(linear)

# --- 1. System Parameters ---
# Transmitter
tx_power_dbm = 30.0
tx_ant_gain_db = 20.0
tx_ant_loss_db = 1.0
tx_filter_loss_db = 1.0
tx_cable_loss_db = 1.0

# Path
frequency_hz = 24e9  # 24 GHz
distance_m = 10e3    # 10 km

# Receiver
rx_ant_gain_db = 1.0
rx_ant_loss_db = 0.5
signal_bw_hz = 100e3 # 100 kHz
temp_k = 300.0

# Receiver Chain components [Name, Gain (dB), Noise Figure (dB)]
# Note: For passive components, loss (dB) is the same as the noise figure (dB) and gain is -loss.
rx_chain = [
    {'name': 'Rx Input Filter', 'gain_db': -1.0, 'nf_db': 1.0},
    {'name': 'LNA',             'gain_db': 36.0, 'nf_db': 2.0},
    {'name': 'Mixer',           'gain_db': -9.0, 'nf_db': 9.0},
    {'name': 'IF Filter',       'gain_db': -1.0, 'nf_db': 1.0},
    {'name': 'IF Amplifier',    'gain_db': 23.0, 'nf_db': 0.0}, # Negligible NF
    {'name': 'Output Filter',   'gain_db': -1.0, 'nf_db': 1.0},
]

# --- 2. Calculations ---

# Step 1: Calculate Transmitter EIRP
print("Step 1: Calculate Transmitter EIRP")
total_tx_loss_db = tx_ant_loss_db + tx_filter_loss_db + tx_cable_loss_db
eirp_dbm = tx_power_dbm + tx_ant_gain_db - total_tx_loss_db
print(f"EIRP (dBm) = Tx Power ({tx_power_dbm} dBm) + Tx Ant Gain ({tx_ant_gain_db} dB) - Tx Losses ({total_tx_loss_db} dB)")
print(f"EIRP = {eirp_dbm:.2f} dBm\n")


# Step 2: Calculate Free Space Path Loss (FSPL)
print("Step 2: Calculate Free Space Path Loss (FSPL)")
# FSPL (dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
fspl_db = 20 * math.log10(distance_m) + 20 * math.log10(frequency_hz) + 20 * math.log10(4 * math.pi / C_LIGHT)
print(f"FSPL (dB) for {distance_m/1000} km at {frequency_hz/1e9} GHz is {fspl_db:.2f} dB\n")

# Step 3: Calculate Received Signal Power (S) at receiver input
print("Step 3: Calculate Received Signal Power (S)")
# This is the power at the input of the first electronic component (the Rx Input Filter)
received_signal_dbm = eirp_dbm - fspl_db + rx_ant_gain_db - rx_ant_loss_db
print(f"S (dBm) = EIRP ({eirp_dbm:.2f} dBm) - FSPL ({fspl_db:.2f} dB) + Rx Ant Gain ({rx_ant_gain_db} dB) - Rx Ant Loss ({rx_ant_loss_db} dB)")
print(f"S = {received_signal_dbm:.2f} dBm\n")


# Step 4: Calculate Receiver System Noise Figure (NF_sys) using Friis Formula
print("Step 4: Calculate System Noise Figure (NF_sys)")
cumulative_gain_linear = 1.0
nf_sys_linear = 0.0
for i, component in enumerate(rx_chain):
    gain_linear = db_to_linear(component['gain_db'])
    nf_linear = db_to_linear(component['nf_db'])

    if i == 0:
        nf_sys_linear += nf_linear
    else:
        nf_sys_linear += (nf_linear - 1) / cumulative_gain_linear
    
    cumulative_gain_linear *= gain_linear

nf_sys_db = linear_to_db(nf_sys_linear)
print(f"The cascaded Noise Figure of the receiver chain is {nf_sys_db:.2f} dB\n")


# Step 5: Calculate Total Input-Referred Noise Power (N)
print("Step 5: Calculate Total Input-Referred Noise Power (N)")
# Thermal noise N = k * T * B
thermal_noise_watts = K_BOLTZMANN * temp_k * signal_bw_hz
thermal_noise_dbm = linear_to_db(thermal_noise_watts) + 30
print(f"Thermal Noise Floor (dBm) for {signal_bw_hz/1000} kHz BW at {temp_k}K = {thermal_noise_dbm:.2f} dBm")

total_noise_dbm = thermal_noise_dbm + nf_sys_db
print(f"N (dBm) = Thermal Noise ({thermal_noise_dbm:.2f} dBm) + System NF ({nf_sys_db:.2f} dB)")
print(f"N = {total_noise_dbm:.2f} dBm\n")

# Step 6: Calculate Final SNR
print("--- Final SNR Calculation ---")
snr_db = received_signal_dbm - total_noise_dbm
print(f"SNR (dB) = Received Signal Power (S) - Total Input-Referred Noise (N)")
print(f"SNR (dB) = {received_signal_dbm:.2f} dBm - ({total_noise_dbm:.2f} dBm)")
print(f"SNR (dB) = {snr_db:.2f} dB")

print("\n<<<" + f"{snr_db:.2f}" + ">>>")