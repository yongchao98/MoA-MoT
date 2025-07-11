import math

def to_linear(db_value):
    """Converts a dB value to a linear ratio."""
    return 10**(db_value / 10)

def to_db(linear_value):
    """Converts a linear ratio to dB."""
    return 10 * math.log10(linear_value)

# --- Input Parameters ---
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
rx_filter_loss_db = 1.0
lna_gain_db = 36.0
lna_nf_db = 2.0
mixer_loss_db = 9.0
mixer_nf_db = 9.0  # Assume NF is equal to conversion loss
if_filter_loss_db = 1.0
if_amp_gain_db = 23.0
if_amp_nf_db = 0.0 # Negligible
if_output_filter_loss_db = 1.0


# System
signal_bw_hz = 100e3  # 100 kHz
temp_k = 300.0       # 300K

# --- Constants ---
k = 1.380649e-23  # Boltzmann's constant
c = 299792458     # Speed of light

# --- Step 1: Calculate Transmitter EIRP ---
print("--- 1. Transmitter Effective Isotropic Radiated Power (EIRP) ---")
total_tx_loss_db = tx_ant_loss_db + tx_filter_loss_db + tx_cable_loss_db
eirp_dbm = tx_power_dbm + tx_ant_gain_db - total_tx_loss_db
print(f"Tx Power: {tx_power_dbm} dBm")
print(f"Tx Antenna Gain: {tx_ant_gain_db} dB")
print(f"Total Tx Losses (Antenna+Filter+Cable): {total_tx_loss_db:.1f} dB")
print(f"EIRP = {tx_power_dbm} + {tx_ant_gain_db} - {total_tx_loss_db:.1f} = {eirp_dbm:.2f} dBm\n")

# --- Step 2: Calculate Free Space Path Loss (FSPL) ---
print("--- 2. Free Space Path Loss (FSPL) ---")
fspl_db = 20 * math.log10(distance_m) + 20 * math.log10(frequency_hz) + 20 * math.log10(4 * math.pi / c)
print(f"Frequency: {frequency_hz / 1e9} GHz")
print(f"Distance: {distance_m / 1e3} km")
print(f"FSPL = {fspl_db:.2f} dB\n")

# --- Step 3: Calculate Received Signal Power (S) ---
print("--- 3. Received Signal Power (S) at Receiver Input ---")
# This is the power available from the antenna, before any Rx components.
s_dbm = eirp_dbm - fspl_db + rx_ant_gain_db - rx_ant_loss_db
print(f"Signal Power = EIRP - FSPL + Rx Ant Gain - Rx Ant Loss")
print(f"S = {eirp_dbm:.2f} dBm - {fspl_db:.2f} dB + {rx_ant_gain_db} dB - {rx_ant_loss_db} dB = {s_dbm:.2f} dBm\n")

# --- Step 4: Calculate Receiver Total Noise Figure (NF) ---
print("--- 4. Receiver Cascaded Noise Figure (NF_rx) ---")
# The receiver chain for NF calculation starts with the first component after the antenna.
# The NF is referred to the input of this first component.
# Stage 1: Rx Input Filter
g1_lin = to_linear(-rx_filter_loss_db)
f1_lin = to_linear(rx_filter_loss_db) # NF of a passive component is its loss

# Stage 2: LNA
g2_lin = to_linear(lna_gain_db)
f2_lin = to_linear(lna_nf_db)

# Stage 3: Mixer
g3_lin = to_linear(-mixer_loss_db)
f3_lin = to_linear(mixer_nf_db)

# Stage 4: IF Filter
g4_lin = to_linear(-if_filter_loss_db)
f4_lin = to_linear(if_filter_loss_db)

# Friis formula for cascaded noise figure
# F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
f_total_lin = f1_lin + (f2_lin - 1) / g1_lin + (f3_lin - 1) / (g1_lin * g2_lin) + (f4_lin - 1) / (g1_lin * g2_lin * g3_lin)
nf_total_db = to_db(f_total_lin)
print(f"Using Friis formula for the receiver chain (Filter -> LNA -> Mixer -> ...)")
print(f"Total Receiver Noise Figure (NF_rx) = {nf_total_db:.2f} dB\n")


# --- Step 5: Calculate Total Input-Referred Noise Power (N) ---
print("--- 5. Total Input-Referred Noise Power (N) ---")
# a) Calculate thermal noise floor
n_thermal_watts = k * temp_k * signal_bw_hz
n_thermal_dbm = 10 * math.log10(n_thermal_watts / 0.001)
print(f"Thermal Noise (kTB) in {signal_bw_hz/1e3} kHz BW at {temp_k}K = {n_thermal_dbm:.2f} dBm")

# b) Add receiver noise figure
n_total_dbm = n_thermal_dbm + nf_total_db
print(f"Total Noise Power = Thermal Noise + Receiver NF")
print(f"N = {n_thermal_dbm:.2f} dBm + {nf_total_db:.2f} dB = {n_total_dbm:.2f} dBm\n")

# --- Step 6: Calculate Final SNR ---
print("--- 6. Final Signal-to-Noise Ratio (SNR) ---")
snr_db = s_dbm - n_total_dbm
print(f"SNR (dB) = Signal Power (dBm) - Total Noise Power (dBm)")
print(f"SNR (dB) = {s_dbm:.2f} - ({n_total_dbm:.2f})")
print(f"Resulting SNR = {snr_db:.2f} dB")

# --- Final Answer ---
final_answer = round(snr_db, 2)
print(f"\n<<<SNR = {final_answer} dB>>>")