import math

# Define the known values from the problem statement.
signal_level = 75.0  # dB

# Data for the train noise
train_level_ref = 100.0  # dB
train_dist_ref = 10.0    # meters
train_dist_target = 30.0 # meters (distance to the people)

# Data for the construction noise
const_level_ref = 115.0  # dB
const_dist_ref = 20.0    # meters
const_dist_target = 50.0 # meters (distance to the people)

# --- Step 1: Calculate the noise level from the train at the people's location ---
print("Step 1: Calculate the train's noise level at the location.")
# Use the formula: L2 = L1 - 20 * log10(r2 / r1)
train_level_at_location = train_level_ref - 20 * math.log10(train_dist_target / train_dist_ref)
print(f"L_train = {train_level_ref} dB - 20 * log10({train_dist_target:.1f} m / {train_dist_ref:.1f} m) = {train_level_at_location:.2f} dB\n")

# --- Step 2: Calculate the noise level from the construction at the people's location ---
print("Step 2: Calculate the construction's noise level at the location.")
# Use the formula: L2 = L1 - 20 * log10(r2 / r1)
const_level_at_location = const_level_ref - 20 * math.log10(const_dist_target / const_dist_ref)
print(f"L_construction = {const_level_ref} dB - 20 * log10({const_dist_target:.1f} m / {const_dist_ref:.1f} m) = {const_level_at_location:.2f} dB\n")

# --- Step 3: Calculate the total noise level by combining the two sources ---
print("Step 3: Calculate the total noise level by summing the intensities.")
# Use the formula: L_total = 10 * log10(10^(L_train/10) + 10^(L_const/10))
total_noise_level = 10 * math.log10(10**(train_level_at_location / 10) + 10**(const_level_at_location / 10))
print(f"L_total_noise = 10 * log10(10^({train_level_at_location:.2f}/10) + 10^({const_level_at_location:.2f}/10)) = {total_noise_level:.2f} dB\n")

# --- Step 4: Calculate the final Signal-to-Noise Ratio (SNR) ---
print("Step 4: Calculate the Signal-to-Noise Ratio (SNR).")
# Use the formula: SNR = Signal - Noise
snr = signal_level - total_noise_level
print(f"SNR = {signal_level:.1f} dB - {total_noise_level:.2f} dB = {snr:.2f} dB")