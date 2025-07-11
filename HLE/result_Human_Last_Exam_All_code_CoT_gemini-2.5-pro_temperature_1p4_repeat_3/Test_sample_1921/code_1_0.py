import math

# --- Given Parameters ---
# Transmitter
P_tx_dbm = 30.0         # Transmitter power in dBm
G_tx_ant_db = 20.0      # Tx antenna gain in dB
L_tx_ant_db = 1.0       # Tx antenna loss in dB
L_tx_filt_db = 1.0      # Tx filter loss in dB
L_tx_cable_db = 1.0     # Tx cable loss in dB
f_hz = 24e9             # Frequency in Hz (24 GHz)
B_hz = 100e3            # Signal bandwidth in Hz (100 kHz)

# Path
dist_m = 10e3           # Distance in meters (10 km)

# Receiver
G_rx_ant_db = 1.0       # Rx antenna gain in dB
L_rx_ant_db = 0.5       # Rx antenna loss in dB
L_rx_filt1_db = 1.0     # Rx input filter loss in dB
G_lna_db = 36.0         # LNA gain in dB
NF_lna_db = 2.0         # LNA Noise Figure in dB
L_mixer_db = 9.0        # Mixer conversion loss in dB (and its NF)
L_rx_filt2_db = 1.0     # IF filter loss in dB
# Negligible NF for subsequent IF amplifier means it doesn't contribute to noise figure

# Constants
T_k = 300.0             # Ambient temperature in Kelvin
k_boltzmann = 1.380649e-23  # Boltzmann's constant in J/K
c_light = 299792458.0    # Speed of light in m/s

# --- Helper functions ---
def db_to_linear(db_val):
    return 10**(db_val / 10.0)

def linear_to_db(lin_val):
    return 10 * math.log10(lin_val)

# --- Step-by-step Calculation ---

print("### Step 1: Calculate Transmitter EIRP ###")
tx_losses_db = L_tx_filt_db + L_tx_cable_db
power_to_ant_dbm = P_tx_dbm - tx_losses_db
eirp_dbm = power_to_ant_dbm + G_tx_ant_db - L_tx_ant_db
print("EIRP (dBm) = Tx Power (dBm) - Tx Losses (dB) + Tx Antenna Gain (dB) - Tx Antenna Loss (dB)")
print(f"EIRP (dBm) = {P_tx_dbm} - ({L_tx_filt_db} + {L_tx_cable_db}) + {G_tx_ant_db} - {L_tx_ant_db} = {eirp_dbm:.2f} dBm\n")


print("### Step 2: Calculate Free Space Path Loss (FSPL) ###")
# FSPL (dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
fspl_db = 20 * math.log10(dist_m) + 20 * math.log10(f_hz) + 20 * math.log10(4 * math.pi / c_light)
print("FSPL (dB) = 20*log10(distance_m) + 20*log10(frequency_hz) + 20*log10(4*pi/c)")
print(f"FSPL (dB) = 20*log10({dist_m:.0f}) + 20*log10({f_hz:.0e}) + 20*log10(4*pi/{c_light:.0e}) = {fspl_db:.2f} dB\n")


print("### Step 3: Calculate Received Signal Power (S_in) ###")
# This is the power at the input of the receiver chain (after Rx antenna)
p_rx_iso_dbm = eirp_dbm - fspl_db
s_in_dbm = p_rx_iso_dbm + G_rx_ant_db - L_rx_ant_db
print("S_in (dBm) = EIRP (dBm) - FSPL (dB) + Rx Antenna Gain (dB) - Rx Antenna Loss (dB)")
print(f"S_in (dBm) = {eirp_dbm:.2f} - {fspl_db:.2f} + {G_rx_ant_db} - {L_rx_ant_db} = {s_in_dbm:.2f} dBm\n")


print("### Step 4: Calculate System Noise Figure (NF_sys) ###")
# We use the Friis formula for cascaded noise figure.
# F_sys = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
# For passive components, NF = Loss.
# Components: 1=Input Filter, 2=LNA, 3=Mixer, 4=IF Filter
G1_lin = db_to_linear(-L_rx_filt1_db)
F1_lin = db_to_linear(L_rx_filt1_db)
G2_lin = db_to_linear(G_lna_db)
F2_lin = db_to_linear(NF_lna_db)
G3_lin = db_to_linear(-L_mixer_db)
F3_lin = db_to_linear(L_mixer_db)
G4_lin = db_to_linear(-L_rx_filt2_db)
F4_lin = db_to_linear(L_rx_filt2_db)

F_sys_lin = F1_lin + (F2_lin - 1)/G1_lin + (F3_lin - 1)/(G1_lin*G2_lin) + (F4_lin-1)/(G1_lin*G2_lin*G3_lin)
nf_sys_db = linear_to_db(F_sys_lin)
print("The system noise figure is calculated by cascading the contributions of each component in the receiver chain.")
print("NF_sys (dB) = 10*log10(F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...)")
print(f"Calculated System Noise Figure (NF_sys) = {nf_sys_db:.2f} dB\n")


print("### Step 5: Calculate Input-Referred Noise Power (N_in) ###")
# N_in (dBm) = 10*log10(k * T * 1000) + 10*log10(B) + NF_sys
n_thermal_density_dbm_hz = 10 * math.log10(k_boltzmann * T_k * 1000)
n_in_dbm = n_thermal_density_dbm_hz + 10 * math.log10(B_hz) + nf_sys_db
print("N_in (dBm) = Thermal Noise Density (dBm/Hz) + 10*log10(Bandwidth_Hz) + NF_sys (dB)")
print(f"N_in (dBm) = {n_thermal_density_dbm_hz:.2f} + {10*math.log10(B_hz):.2f} + {nf_sys_db:.2f} = {n_in_dbm:.2f} dBm\n")


print("### Step 6: Calculate Final Signal-to-Noise Ratio (SNR) ###")
snr_db = s_in_dbm - n_in_dbm
print("Final SNR (dB) = S_in (dBm) - N_in (dBm)")
print(f"Final SNR (dB) = {s_in_dbm:.2f} dBm - ({n_in_dbm:.2f} dBm) = {snr_db:.2f} dB")
print("<<<28.27>>>")