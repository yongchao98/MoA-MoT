import math

# Function to convert dB to linear scale for gain
def db_to_linear(db_value):
    return 10**(db_value / 10)

# Function to convert linear scale to dB
def linear_to_db(linear_value):
    return 10 * math.log10(linear_value)

# 0. Define Constants
tx_power_dbm = 30.0
tx_antenna_gain_db = 20.0
tx_antenna_loss_db = 1.0
tx_filter_loss_db = 1.0
tx_cable_loss_db = 1.0
frequency_hz = 24e9
distance_m = 10e3
rx_antenna_gain_db = 1.0
rx_antenna_loss_db = 0.5
rx_input_filter_loss_db = 1.0
lna_gain_db = 36.0
lna_nf_db = 2.0
mixer_loss_db = 9.0
if_filter_loss_db = 1.0
if_amp_gain_db = 23.0
output_filter_loss_db = 1.0
temperature_k = 300.0
bandwidth_hz = 100e3
k_boltzmann = 1.38e-23

# 1. Calculate Effective Isotropic Radiated Power (EIRP)
total_tx_loss_db = tx_antenna_loss_db + tx_filter_loss_db + tx_cable_loss_db
eirp_dbm = tx_power_dbm + tx_antenna_gain_db - total_tx_loss_db
print(f"Step 1: Transmitter EIRP Calculation")
print(f"  EIRP (dBm) = Tx Power (dBm) + Tx Gain (dB) - Tx Losses (dB)")
print(f"  EIRP (dBm) = {tx_power_dbm} + {tx_antenna_gain_db} - ({tx_antenna_loss_db} + {tx_filter_loss_db} + {tx_cable_loss_db}) = {eirp_dbm:.2f} dBm\n")


# 2. Calculate Free Space Path Loss (FSPL)
# FSPL (dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
# 20*log10(4*pi/c) is approximately -147.55 for frequency in Hz and distance in m.
fspl_db = 20 * math.log10(distance_m) + 20 * math.log10(frequency_hz) - 147.55
print(f"Step 2: Free Space Path Loss Calculation")
print(f"  FSPL for {distance_m/1000} km at {frequency_hz/1e9} GHz = {fspl_db:.2f} dB\n")

# 3. Calculate Received Signal Power (S_in) at the input reference plane
# The reference plane is before the Rx antenna's own loss.
rx_power_at_antenna_dbm = eirp_dbm - fspl_db + rx_antenna_gain_db
print(f"Step 3: Received Signal Power (S) Calculation")
print(f"  S (dBm) = EIRP (dBm) - FSPL (dB) + Rx Gain (dB)")
print(f"  S (dBm) = {eirp_dbm:.2f} - {fspl_db:.2f} + {rx_antenna_gain_db} = {rx_power_at_antenna_dbm:.2f} dBm\n")

# 4. Calculate Receiver System Noise Figure (NF_sys)
# Loss before LNA adds directly to the noise figure in dB.
loss_before_lna_db = rx_antenna_loss_db + rx_input_filter_loss_db

# Use Friis formula for the cascade starting from the LNA
# F_cascade = F_lna + (F_mixer-1)/G_lna + (F_if_filter-1)/(G_lna*G_mixer) + ...
g_lna = db_to_linear(lna_gain_db)
f_lna = db_to_linear(lna_nf_db)
g_mixer = db_to_linear(-mixer_loss_db) # Loss is negative gain
f_mixer = db_to_linear(mixer_loss_db) # For passive mixers, NF = Loss
f_if_filter = db_to_linear(if_filter_loss_db)

f_cascade = f_lna + (f_mixer - 1) / g_lna + (f_if_filter - 1) / (g_lna * g_mixer)
nf_cascade_db = linear_to_db(f_cascade)

# Total system NF is the sum of pre-LNA loss and the subsequent cascade NF
nf_sys_db = loss_before_lna_db + nf_cascade_db
print(f"Step 4: Receiver System Noise Figure (NF) Calculation")
print(f"  NF_sys (dB) = Loss_pre_LNA (dB) + NF_cascade (dB)")
print(f"  NF_sys (dB) = ({rx_antenna_loss_db} + {rx_input_filter_loss_db}) + {nf_cascade_db:.2f} = {nf_sys_db:.2f} dB\n")

# 5. Calculate Total Input-Referred Noise Power (N_in)
# Calculate thermal noise in the given bandwidth
noise_thermal_watts = k_boltzmann * temperature_k * bandwidth_hz
noise_thermal_dbm = 10 * math.log10(noise_thermal_watts / 0.001)

# Input-referred noise is thermal noise + system noise figure
total_input_noise_dbm = noise_thermal_dbm + nf_sys_db
print(f"Step 5: Input-Referred Noise Power (N) Calculation")
print(f"  N (dBm) = Thermal Noise (dBm) + NF_sys (dB)")
print(f"  N (dBm) = {noise_thermal_dbm:.2f} + {nf_sys_db:.2f} = {total_input_noise_dbm:.2f} dBm\n")


# 6. Calculate Final SNR
snr_db = rx_power_at_antenna_dbm - total_input_noise_dbm
print(f"Step 6: Final Signal-to-Noise Ratio (SNR) Calculation")
print(f"The final equation for SNR is: Signal Power (dBm) - Noise Power (dBm)")
print(f"Resulting SNR = {rx_power_at_antenna_dbm:.2f} dBm - ({total_input_noise_dbm:.2f} dBm)")
print(f"Final Result = {snr_db:.2f} dB")
print(f"<<<{snr_db:.2f}>>>")