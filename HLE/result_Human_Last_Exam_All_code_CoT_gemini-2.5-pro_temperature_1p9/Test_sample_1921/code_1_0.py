import math

def db_to_linear(db_value):
    """Converts a value in dB to a linear scale."""
    return 10**(db_value / 10)

def linear_to_db(linear_value):
    """Converts a linear value to dB."""
    return 10 * math.log10(linear_value)

# 1. System Parameters
# Constants
k = 1.38e-23  # Boltzmann's constant in J/K
T = 300       # Ambient temperature in Kelvin
B = 100e3     # Signal bandwidth in Hz

# Transmitter (Tx) Parameters
P_tx_dBm = 30.0
G_tx_ant_dB = 20.0
L_tx_ant_dB = 1.0
L_tx_filter_dB = 1.0
L_tx_cable_dB = 1.0

# Path Parameters
f_GHz = 24.0
d_km = 10.0

# Receiver (Rx) Parameters
G_rx_ant_dB = 1.0
L_rx_ant_dB = 0.5
L_rx_filter_dB = 1.0

# Receiver Chain for Noise Figure Calculation (starting from LNA input)
lna_gain_dB = 36.0
lna_nf_dB = 2.0
mixer_loss_dB = 9.0
if_filter_loss_dB = 1.0
if_amp_gain_dB = 23.0
if_amp_nf_dB = 0.0 # Negligible
out_filter_loss_dB = 1.0

print("Calculating Signal-to-Noise Ratio (SNR)\n")

# 2. Calculate Effective Isotropic Radiated Power (EIRP)
print("--- Step 1: Calculate Transmitter EIRP ---")
L_tx_total_dB = L_tx_ant_dB + L_tx_filter_dB + L_tx_cable_dB
EIRP_dBm = P_tx_dBm + G_tx_ant_dB - L_tx_total_dB
print(f"EIRP = Tx Power [{P_tx_dBm:.1f} dBm] + Tx Gain [{G_tx_ant_dB:.1f} dB] - Tx Losses [{L_tx_total_dB:.1f} dB] = {EIRP_dBm:.2f} dBm\n")

# 3. Calculate Free Space Path Loss (FSPL)
print("--- Step 2: Calculate Free Space Path Loss (FSPL) ---")
FSPL_dB = 20 * math.log10(d_km) + 20 * math.log10(f_GHz) + 92.45
print(f"FSPL for {d_km:.1f} km at {f_GHz:.1f} GHz = {FSPL_dB:.2f} dB\n")

# 4. Calculate Received Signal Power at LNA input (Prx)
print("--- Step 3: Calculate Received Signal Power (Prx) at LNA Input ---")
L_pre_LNA_dB = L_rx_ant_dB + L_rx_filter_dB
P_rx_at_LNA_dBm = EIRP_dBm - FSPL_dB + G_rx_ant_dB - L_pre_LNA_dB
print(f"Prx = EIRP [{EIRP_dBm:.2f} dBm] - FSPL [{FSPL_dB:.2f} dB] + Rx Gain [{G_rx_ant_dB:.1f} dB] - Pre-LNA Losses [{L_pre_LNA_dB:.1f} dB]")
print(f"Signal Power (Prx) = {P_rx_at_LNA_dBm:.2f} dBm\n")

# 5. Calculate Total Noise Power at LNA input
print("--- Step 4: Calculate Total Noise Power (N_total) at LNA Input ---")

# 5a. Calculate Receiver Noise Figure (NF_rx) using Friis formula
# Convert component values to linear scale
# Stage 1: LNA
F1_lin = db_to_linear(lna_nf_dB)
G1_lin = db_to_linear(lna_gain_dB)
# Stage 2: Mixer (NF = Loss, Gain = 1/Loss)
F2_lin = db_to_linear(mixer_loss_dB)
G2_lin = db_to_linear(-mixer_loss_dB)
# Stage 3: IF Filter (NF = Loss, Gain = 1/Loss)
F3_lin = db_to_linear(if_filter_loss_dB)
G3_lin = db_to_linear(-if_filter_loss_dB)
# Stage 4: IF Amplifier
F4_lin = db_to_linear(if_amp_nf_dB) # NF=0dB -> F=1
G4_lin = db_to_linear(if_amp_gain_dB)
# Stage 5: Output Filter (NF = Loss, Gain = 1/Loss)
F5_lin = db_to_linear(out_filter_loss_dB)

# Friis Formula for Noise Factor: F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
F_rx_lin = F1_lin + (F2_lin - 1)/G1_lin + (F3_lin - 1)/(G1_lin*G2_lin) + (F4_lin - 1)/(G1_lin*G2_lin*G3_lin) + (F5_lin - 1)/(G1_lin*G2_lin*G3_lin*G4_lin)
NF_rx_dB = linear_to_db(F_rx_lin)
print(f"Cascaded Receiver Noise Figure (NF_rx) = {NF_rx_dB:.2f} dB")

# 5b. Calculate total noise power
N_thermal_watts = k * T * B
N_thermal_dBm = 10 * math.log10(N_thermal_watts / 0.001)
N_total_dBm = N_thermal_dBm + NF_rx_dB
print(f"Thermal Noise (kTB) in {B/1000:.0f} kHz bandwidth = {N_thermal_dBm:.2f} dBm")
print(f"Total Noise Power (N_total) = Thermal Noise [{N_thermal_dBm:.2f} dBm] + NF_rx [{NF_rx_dB:.2f} dB] = {N_total_dBm:.2f} dBm\n")

# 6. Final SNR Calculation
print("--- Step 5: Final SNR Calculation ---")
SNR_dB = P_rx_at_LNA_dBm - N_total_dBm
print(f"SNR (dB) = Prx (dBm) - N_total (dBm)")
print(f"SNR (dB) = {P_rx_at_LNA_dBm:.2f} dBm - ({N_total_dBm:.2f} dBm)")
print(f"Resulting SNR = {SNR_dB:.2f} dB")

final_answer = f"{SNR_dB:.2f}"
print(f"\n<<<{final_answer}>>>")
