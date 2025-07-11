import math

def db_to_linear(db):
    """Converts a value in dB to a linear ratio."""
    return 10**(db / 10)

def linear_to_db(linear):
    """Converts a linear ratio to a value in dB."""
    if linear <= 0:
        return -float('inf')
    return 10 * math.log10(linear)

# --- System Parameters ---
# Transmitter (Tx)
TX_POWER_DBM = 30.0
TX_ANT_GAIN_DB = 20.0
TX_ANT_LOSS_DB = 1.0
TX_FILTER_LOSS_DB = 1.0
TX_CABLE_LOSS_DB = 1.0

# Path
FREQUENCY_GHZ = 24.0
DISTANCE_KM = 10.0

# Receiver (Rx)
RX_ANT_GAIN_DB = 1.0
RX_ANT_LOSS_DB = 0.5
RX_INPUT_FILTER_LOSS_DB = 1.0
LNA_GAIN_DB = 36.0
LNA_NF_DB = 2.0
MIXER_LOSS_DB = 9.0
MIXER_NF_DB = 9.0  # Noise figure of a passive mixer is its conversion loss
IF_FILTER_LOSS_DB = 1.0
IF_AMP_GAIN_DB = 23.0
IF_AMP_NF_DB = 0.0   # Negligible NF, 0 dB = 1 linear
OUTPUT_FILTER_LOSS_DB = 1.0

# System
BANDWIDTH_KHZ = 100.0
TEMPERATURE_K = 300.0
BOLTZMANN_K = 1.380649e-23

print("--- Step 1: Calculate Effective Isotropic Radiated Power (EIRP) ---")
tx_losses = TX_ANT_LOSS_DB + TX_FILTER_LOSS_DB + TX_CABLE_LOSS_DB
eirp_dbm = TX_POWER_DBM + TX_ANT_GAIN_DB - tx_losses
print(f"EIRP (dBm) = Tx Power (dBm) + Tx Antenna Gain (dB) - Tx Losses (dB)")
print(f"EIRP = {TX_POWER_DBM} + {TX_ANT_GAIN_DB} - ({TX_ANT_LOSS_DB} + {TX_FILTER_LOSS_DB} + {TX_CABLE_LOSS_DB}) = {eirp_dbm:.2f} dBm\n")

print("--- Step 2: Calculate Free Space Path Loss (FSPL) ---")
# FSPL (dB) = 92.45 + 20*log10(f) + 20*log10(d) where f is in GHz and d is in km
fspl_db = 92.45 + 20 * math.log10(FREQUENCY_GHZ) + 20 * math.log10(DISTANCE_KM)
print(f"FSPL (dB) = 92.45 + 20*log10(Frequency_GHz) + 20*log10(Distance_km)")
print(f"FSPL = 92.45 + 20*log10({FREQUENCY_GHZ}) + 20*log10({DISTANCE_KM}) = {fspl_db:.2f} dB\n")

print("--- Step 3: Calculate Received Signal Power (S) ---")
# Signal power at the input of the receiver cascade (after Rx antenna, before Rx filter)
rx_power_dbm = eirp_dbm - fspl_db + RX_ANT_GAIN_DB - RX_ANT_LOSS_DB
print(f"Signal Power (dBm) = EIRP (dBm) - FSPL (dB) + Rx Antenna Gain (dB) - Rx Antenna Loss (dB)")
print(f"S = {eirp_dbm:.2f} - {fspl_db:.2f} + {RX_ANT_GAIN_DB} - {RX_ANT_LOSS_DB} = {rx_power_dbm:.2f} dBm\n")

print("--- Step 4: Calculate Receiver Cascaded Noise Figure (NF_rx) ---")
# Cascade starts at the input of the Rx input filter. Friis formula: F_total = F1 + (F2-1)/G1 + ...
stages = [
    {'name': 'Rx Input Filter', 'gain_db': -RX_INPUT_FILTER_LOSS_DB, 'nf_db': RX_INPUT_FILTER_LOSS_DB},
    {'name': 'LNA',             'gain_db': LNA_GAIN_DB,              'nf_db': LNA_NF_DB},
    {'name': 'Mixer',           'gain_db': -MIXER_LOSS_DB,           'nf_db': MIXER_NF_DB},
    {'name': 'IF Filter',       'gain_db': -IF_FILTER_LOSS_DB,       'nf_db': IF_FILTER_LOSS_DB},
    {'name': 'IF Amplifier',    'gain_db': IF_AMP_GAIN_DB,           'nf_db': IF_AMP_NF_DB},
    {'name': 'Output Filter',   'gain_db': -OUTPUT_FILTER_LOSS_DB,   'nf_db': OUTPUT_FILTER_LOSS_DB}
]

total_nf_linear = 0
cumulative_gain_linear = 1
print("Receiver NF calculated using Friis formula. Contribution of each stage:")
for i, stage in enumerate(stages):
    gain_linear = db_to_linear(stage['gain_db'])
    nf_linear = db_to_linear(stage['nf_db'])
    stage_contribution = (nf_linear - 1) / cumulative_gain_linear
    if i == 0: # First stage contribution is F1-1, not (F1-1)/G_previous
        stage_contribution = nf_linear 
        total_nf_linear += stage_contribution - 1 # Adjust for F_total=F1+(F2-1)/G1.. formula
    
    total_nf_linear += stage_contribution
    cumulative_gain_linear *= gain_linear
    print(f"Term {i+1} ({stage['name']}): {stage_contribution:.5f}")

rx_nf_db = linear_to_db(total_nf_linear)
print(f"\nTotal Receiver Noise Figure (Linear) = {total_nf_linear:.4f}")
print(f"Total Receiver Noise Figure (dB) = {rx_nf_db:.2f} dB\n")

print("--- Step 5: Calculate Total Noise Power (N) ---")
bandwidth_hz = BANDWIDTH_KHZ * 1000
# Thermal noise density at T Kelvin, in dBm/Hz
noise_density_dbm_hz = linear_to_db(BOLTZMANN_K * TEMPERATURE_K * 1000)
# Total Noise Power in dBm
noise_power_dbm = noise_density_dbm_hz + 10 * math.log10(bandwidth_hz) + rx_nf_db
print(f"Noise Power (dBm) = Thermal Noise (dBm/Hz) + 10*log10(Bandwidth_Hz) + NF_rx (dB)")
print(f"N = {noise_density_dbm_hz:.2f} dBm/Hz + {10 * math.log10(bandwidth_hz):.2f} dB + {rx_nf_db:.2f} dB = {noise_power_dbm:.2f} dBm\n")


print("--- Step 6: Calculate Final Signal-to-Noise Ratio (SNR) ---")
snr_db = rx_power_dbm - noise_power_dbm
print(f"SNR (dB) = Signal Power (dBm) - Noise Power (dBm)")
print(f"SNR = {rx_power_dbm:.2f} dBm - ({noise_power_dbm:.2f} dBm) = {snr_db:.2f} dB")

print(f"\nFinal Resulting SNR: {snr_db:.2f} dB")
<<<28.27>>>