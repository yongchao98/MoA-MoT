import math

# Helper functions to convert between dB and linear scale
def db_to_linear(db):
    """Converts a dB value to a linear ratio."""
    return 10**(db / 10)

def linear_to_db(linear):
    """Converts a linear ratio to a dB value."""
    return 10 * math.log10(linear)

# --- Given Parameters ---
# Constants
k = 1.380649e-23  # Boltzmann's constant in J/K
T = 300            # Ambient temperature in Kelvin
B = 100e3          # Signal bandwidth in Hz

# Transmitter (Tx)
tx_power_dbm = 30
tx_ant_gain_db = 20
tx_ant_loss_db = 1
tx_filter_loss_db = 1
tx_cable_loss_db = 1

# Path
freq_ghz = 24
dist_km = 10

# Receiver (Rx) Chain components for Noise Figure calculation
# The chain for NF calculation stops where the cumulative gain is very high
# or when a component has negligible NF (like the IF Amp).
# For passive components (losses), Noise Figure (dB) = Loss (dB).
# Gain (dB) for a lossy component is negative.
rx_chain = [
    {'name': 'Rx Antenna Loss', 'gain_db': -0.5, 'nf_db': 0.5},
    {'name': 'Rx Input Filter', 'gain_db': -1.0, 'nf_db': 1.0},
    {'name': 'LNA',             'gain_db': 36.0, 'nf_db': 2.0},
    {'name': 'Mixer',           'gain_db': -9.0, 'nf_db': 9.0},
    {'name': 'IF Filter',       'gain_db': -1.0, 'nf_db': 1.0},
]
rx_ant_gain_db = 1.0
rx_ant_loss_db = 0.5

# --- Step 1: Calculate Effective Isotropic Radiated Power (EIRP) ---
total_tx_loss_db = tx_ant_loss_db + tx_filter_loss_db + tx_cable_loss_db
eirp_dbm = tx_power_dbm + tx_ant_gain_db - total_tx_loss_db

# --- Step 2: Calculate Free Space Path Loss (FSPL) ---
# FSPL (dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
# Simplified for d in km and f in GHz: FSPL (dB) = 20*log10(d) + 20*log10(f) + 92.45
fspl_db = 20 * math.log10(dist_km) + 20 * math.log10(freq_ghz) + 92.45

# --- Step 3: Calculate Received Signal Power (S_in) ---
# Power at the input to the antenna
power_at_rx_ant_input_dbm = eirp_dbm - fspl_db + rx_ant_gain_db
# Signal power at the reference plane (after antenna loss)
s_in_dbm = power_at_rx_ant_input_dbm - rx_ant_loss_db

# --- Step 4: Calculate System Noise Figure (NF_sys) ---
total_gain_linear = 1.0
total_noise_factor_linear = 0.0

for component in rx_chain:
    gain_linear = db_to_linear(component['gain_db'])
    noise_factor_linear = db_to_linear(component['nf_db'])
    
    # Friis formula for cascaded noise figure
    if total_gain_linear == 1.0: # First component
        component_noise_contribution = noise_factor_linear
    else:
        component_noise_contribution = (noise_factor_linear - 1) / total_gain_linear
        
    total_noise_factor_linear += component_noise_contribution
    total_gain_linear *= gain_linear

nf_sys_db = linear_to_db(total_noise_factor_linear)

# --- Step 5: Calculate Total Input-Referred Noise (N_in) ---
# Thermal noise in watts, then convert to dBm
thermal_noise_watts = k * T * B
thermal_noise_dbm = 10 * math.log10(thermal_noise_watts / 0.001)

# Total input-referred noise is thermal noise + system noise figure
n_in_dbm = thermal_noise_dbm + nf_sys_db

# --- Step 6: Calculate Final SNR ---
snr_db = s_in_dbm - n_in_dbm

# --- Print the results step-by-step ---
print("Calculating SNR of the Received Signal\n")
print(f"1. Transmitter EIRP Calculation:")
print(f"   EIRP (dBm) = Tx Power (dBm) + Tx Ant Gain (dB) - Tx Losses (dB)")
print(f"   EIRP = {tx_power_dbm} dBm + {tx_ant_gain_db} dB - {total_tx_loss_db:.2f} dB = {eirp_dbm:.2f} dBm\n")

print(f"2. Free Space Path Loss Calculation:")
print(f"   FSPL (dB) = 20*log10(Dist_km) + 20*log10(Freq_GHz) + 92.45")
print(f"   FSPL = 20*log10({dist_km}) + 20*log10({freq_ghz}) + 92.45 = {fspl_db:.2f} dB\n")

print(f"3. Received Signal Power (S_in) Calculation:")
print(f"   S_in (dBm) = EIRP (dBm) - FSPL (dB) + Rx Ant Gain (dB) - Rx Ant Loss (dB)")
print(f"   S_in = {eirp_dbm:.2f} dBm - {fspl_db:.2f} dB + {rx_ant_gain_db:.2f} dB - {rx_ant_loss_db:.2f} dB = {s_in_dbm:.2f} dBm\n")

print(f"4. Total Input-Referred Noise (N_in) Calculation:")
print(f"   Thermal Noise Floor (dBm) = 10*log10(k*T*B / 1mW) = {thermal_noise_dbm:.2f} dBm (for {B/1000} kHz BW)")
print(f"   Cascaded System Noise Figure (NF) = {nf_sys_db:.2f} dB")
print(f"   N_in (dBm) = Thermal Noise Floor (dBm) + System NF (dB)")
print(f"   N_in = {thermal_noise_dbm:.2f} dBm + {nf_sys_db:.2f} dB = {n_in_dbm:.2f} dBm\n")

print(f"5. Final SNR Calculation:")
print(f"   SNR (dB) = S_in (dBm) - N_in (dBm)")
print(f"   SNR (dB) = {s_in_dbm:.2f} dBm - ({n_in_dbm:.2f} dBm)")
print(f"   Resulting SNR = {snr_db:.2f} dB")

print(f"\n<<<SNR = {snr_db:.2f} dB>>>")
