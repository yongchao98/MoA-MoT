import math

# --- Helper Functions ---
def db_to_linear(db):
    """Converts a dB value to a linear scale."""
    return 10**(db / 10)

def linear_to_db(linear):
    """Converts a linear value to a dB scale."""
    return 10 * math.log10(linear)

def calculate_fspl(freq_hz, dist_m):
    """Calculates Free Space Path Loss in dB."""
    c = 299792458  # Speed of light in m/s
    # FSPL formula in dB: 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
    fspl_db = 20 * math.log10(dist_m) + 20 * math.log10(freq_hz) + 20 * math.log10((4 * math.pi) / c)
    return fspl_db

# --- Given Parameters ---

# Transmitter
tx_power_dbm = 30.0
tx_antenna_gain_db = 20.0
tx_antenna_loss_db = 1.0
tx_filter_loss_db = 1.0
tx_cable_loss_db = 1.0

# Path
frequency_ghz = 24.0
distance_km = 10.0
signal_bandwidth_khz = 100.0

# Receiver
rx_antenna_gain_db = 1.0
# Receiver components in order (Gain in dB, Noise Figure in dB)
rx_components = [
    {'name': 'Rx Antenna Loss', 'gain_db': -0.5, 'nf_db': 0.5},
    {'name': 'Rx Input Filter', 'gain_db': -1.0, 'nf_db': 1.0},
    {'name': 'LNA',             'gain_db': 36.0, 'nf_db': 2.0},
    {'name': 'Mixer',           'gain_db': -9.0, 'nf_db': 9.0}, # Conversion loss is negative gain, NF equals loss
    {'name': 'IF Filter',       'gain_db': -1.0, 'nf_db': 1.0},
    {'name': 'IF Amplifier',    'gain_db': 23.0, 'nf_db': 0.0}, # Negligible NF
    {'name': 'Output Filter',   'gain_db': -1.0, 'nf_db': 1.0},
]

# System
temperature_k = 300.0

# --- Calculations ---

# 1. Calculate Transmitter EIRP
print("Step 1: Calculating Transmitter EIRP...")
tx_total_loss_db = tx_antenna_loss_db + tx_filter_loss_db + tx_cable_loss_db
eirp_dbm = tx_power_dbm + tx_antenna_gain_db - tx_total_loss_db
print(f"EIRP (dBm) = {tx_power_dbm} dBm + {tx_antenna_gain_db} dB - {tx_total_loss_db} dB = {eirp_dbm:.2f} dBm\n")

# 2. Calculate Free Space Path Loss
print("Step 2: Calculating Free Space Path Loss...")
frequency_hz = frequency_ghz * 1e9
distance_m = distance_km * 1e3
fspl_db = calculate_fspl(frequency_hz, distance_m)
print(f"FSPL for {frequency_ghz} GHz at {distance_km} km = {fspl_db:.2f} dB\n")

# 3. Calculate Received Signal Power (S)
# The reference point for SNR is the input of the receiver, right after the antenna element's gain is applied.
print("Step 3: Calculating Received Signal Power (S)...")
signal_power_dbm = eirp_dbm - fspl_db + rx_antenna_gain_db
print(f"Received Signal Power (dBm) = EIRP - FSPL + Rx Ant Gain")
print(f"S = {eirp_dbm:.2f} dBm - {fspl_db:.2f} dB + {rx_antenna_gain_db} dB = {signal_power_dbm:.2f} dBm\n")

# 4. Calculate Total Receiver Noise (N)
# 4a. Calculate cascaded Noise Figure using Friis Formula
print("Step 4: Calculating Total Receiver Noise Power (N)...")
F_cascade_linear = 0
G_cascade_linear = 1.0
for component in rx_components:
    gain_linear = db_to_linear(component['gain_db'])
    nf_linear = db_to_linear(component['nf_db'])
    
    # Friis formula: F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
    F_cascade_linear += (nf_linear - 1) / G_cascade_linear
    G_cascade_linear *= gain_linear

# Add F1 (the first component's NF) to complete the formula
F_cascade_linear += db_to_linear(rx_components[0]['nf_db'])
nf_cascade_db = linear_to_db(F_cascade_linear)
print(f"Total Cascaded Receiver Noise Figure = {nf_cascade_db:.2f} dB")

# 4b. Calculate Thermal Noise Power
k = 1.380649e-23  # Boltzmann's constant J/K
bandwidth_hz = signal_bandwidth_khz * 1e3
# Thermal Noise in dBm = 10*log10(k*T*B / 1mW)
# A common shortcut: -174 dBm/Hz + 10*log10(Bandwidth) at T=290K. We'll calculate it directly for 300K.
noise_power_watts = k * temperature_k * bandwidth_hz
noise_power_dbm_thermal = 10 * math.log10(noise_power_watts / 0.001)

# 4c. Calculate Total Input Noise Power
total_noise_power_dbm = noise_power_dbm_thermal + nf_cascade_db
print(f"Thermal Noise in {signal_bandwidth_khz} kHz BW = {noise_power_dbm_thermal:.2f} dBm")
print(f"Total Equivalent Input Noise (dBm) = Thermal Noise + Noise Figure")
print(f"N = {noise_power_dbm_thermal:.2f} dBm + {nf_cascade_db:.2f} dB = {total_noise_power_dbm:.2f} dBm\n")

# 5. Calculate Final SNR
print("Step 5: Calculating Final Signal-to-Noise Ratio (SNR)...")
snr_db = signal_power_dbm - total_noise_power_dbm
print("SNR (dB) = Signal Power (dBm) - Total Noise Power (dBm)")
print(f"SNR (dB) = {signal_power_dbm:.2f} dBm - ({total_noise_power_dbm:.2f} dBm) = {snr_db:.2f} dB")

<<<28.46>>>