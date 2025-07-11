import math

def db_to_linear(db):
    """Converts a dB value to a linear scale."""
    return 10**(db / 10)

def linear_to_db(lin):
    """Converts a linear value to a dB scale."""
    return 10 * math.log10(lin)

# Step 1: Define constants and input parameters
# Transmitter (Tx) side
tx_power_dbm = 30.0  # dBm
tx_ant_gain_db = 20.0  # dB
tx_ant_loss_db = 1.0   # dB
tx_filter_loss_db = 1.0 # dB
tx_cable_loss_db = 1.0  # dB

# Path
frequency_ghz = 24.0   # GHz
distance_km = 10.0     # km

# Receiver (Rx) side
rx_ant_gain_db = 1.0   # dB
rx_ant_loss_db = 0.5   # dB
# Receiver components for Noise Figure calculation: (Gain dB, Noise Figure dB)
# Note: For passive components, NF = Loss. Gain = -Loss.
# For negligible NF, NF is 0 dB.
rx_components = [
    (-1.0, 1.0),   # 1. Rx input 150 MHz filter (1 dB loss)
    (36.0, 2.0),   # 2. LNA (36 dB gain, 2 dB NF)
    (-9.0, 9.0),   # 3. Mixer (9 dB conversion loss)
    (-1.0, 1.0),   # 4. 10 MHz IF filter (1 dB loss)
    (23.0, 0.0),   # 5. IF amplifier (23 dB gain, negligible NF)
    (-1.0, 1.0)    # 6. Output filter of 20 MHz (1 dB loss)
]

# System parameters
signal_bw_khz = 100.0  # kHz
temp_k = 300.0         # Kelvin
boltzmann_k = 1.380649e-23  # J/K

# Step 2: Calculate Transmitter EIRP
total_tx_loss_db = tx_ant_loss_db + tx_filter_loss_db + tx_cable_loss_db
eirp_dbm = tx_power_dbm + tx_ant_gain_db - total_tx_loss_db

# Step 3: Calculate Free Space Path Loss (FSPL)
# FSPL (dB) = 32.45 + 20*log10(f_MHz) + 20*log10(d_km)
frequency_mhz = frequency_ghz * 1000
fspl_db = 32.45 + 20 * math.log10(frequency_mhz) + 20 * math.log10(distance_km)

# Step 4: Calculate Received Signal Power (S) at the receiver input
net_rx_ant_gain_db = rx_ant_gain_db - rx_ant_loss_db
signal_power_dbm = eirp_dbm - fspl_db + net_rx_ant_gain_db

# Step 5: Calculate Receiver System's Cascaded Noise Figure (NF_sys)
total_gain_linear = 1.0
total_noise_factor_linear = 0.0

for i, (gain_db, nf_db) in enumerate(rx_components):
    gain_linear = db_to_linear(gain_db)
    noise_factor_linear = db_to_linear(nf_db)
    
    if i == 0:
        total_noise_factor_linear = noise_factor_linear
    else:
        total_noise_factor_linear += (noise_factor_linear - 1) / total_gain_linear
        
    total_gain_linear *= gain_linear

nf_sys_db = linear_to_db(total_noise_factor_linear)

# Step 6: Calculate Total Noise Power (N)
# N0 (dBm/Hz) = 10*log10(k*T*1000)
noise_density_dbm_hz = 10 * math.log10(boltzmann_k * temp_k * 1000)
# Noise in Bandwidth (dB) = 10*log10(BW in Hz)
signal_bw_hz = signal_bw_khz * 1000
noise_in_bw_db = 10 * math.log10(signal_bw_hz)

# Total noise = thermal noise + receiver-generated noise
total_noise_power_dbm = noise_density_dbm_hz + noise_in_bw_db + nf_sys_db

# Step 7: Calculate the final SNR
snr_db = signal_power_dbm - total_noise_power_dbm

# Final Output
print("Calculation Steps:")
print(f"1. Tx EIRP = {tx_power_dbm:.2f}dBm(Tx) + {tx_ant_gain_db:.2f}dB(Gain) - {total_tx_loss_db:.2f}dB(Loss) = {eirp_dbm:.2f} dBm")
print(f"2. FSPL = {fspl_db:.2f} dB (for {frequency_ghz} GHz at {distance_km} km)")
print(f"3. Signal Power (S) = {eirp_dbm:.2f}dBm(EIRP) - {fspl_db:.2f}dB(FSPL) + {net_rx_ant_gain_db:.2f}dB(RxGain) = {signal_power_dbm:.2f} dBm")
print(f"4. Receiver Noise Figure (NF_sys) = {nf_sys_db:.2f} dB")
print(f"5. Noise Power (N) = {noise_density_dbm_hz:.2f}dBm/Hz(N0) + {noise_in_bw_db:.2f}dB(BW) + {nf_sys_db:.2f}dB(NF) = {total_noise_power_dbm:.2f} dBm")
print("\n--- Final Result ---")
print(f"SNR (dB) = S (dBm) - N (dBm)")
print(f"SNR (dB) = {signal_power_dbm:.2f} - ({total_noise_power_dbm:.2f})")
print(f"SNR (dB) = {snr_db:.2f}")
