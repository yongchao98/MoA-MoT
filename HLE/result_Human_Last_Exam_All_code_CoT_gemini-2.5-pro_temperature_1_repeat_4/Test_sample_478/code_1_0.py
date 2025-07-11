import math

# Plan:
# 1. Define the signal level.
# 2. Calculate the noise level from the train at the people's location.
# 3. Calculate the noise level from the construction at the people's location.
# 4. Combine the two noise levels to get the total noise level.
# 5. Calculate the Signal-to-Noise Ratio (SNR).

# --- Given Values ---
# Signal Level at the people's location
L_signal = 75.0  # dB

# Train noise data
L_train_at_10m = 100.0  # dB
d_train_to_people = 30.0  # meters

# Construction noise data
L_const_at_20m = 115.0  # dB
d_const_to_people = 50.0  # meters

# --- Calculations ---

# Step 2: Calculate train noise level at the people's location (30m)
# Formula: L2 = L1 - 20 * log10(r2 / r1)
L_train_at_people = L_train_at_10m - 20 * math.log10(d_train_to_people / 10.0)

# Step 3: Calculate construction noise level at the people's location (50m)
L_const_at_people = L_const_at_20m - 20 * math.log10(d_const_to_people / 20.0)

# Step 4: Combine the two noise levels
# Convert dB to intensity (I = 10^(L/10))
I_train = 10**(L_train_at_people / 10)
I_const = 10**(L_const_at_people / 10)

# Add intensities
I_total_noise = I_train + I_const

# Convert total intensity back to dB (L = 10 * log10(I))
L_total_noise = 10 * math.log10(I_total_noise)

# Step 5: Calculate SNR
# SNR = Signal (dB) - Noise (dB)
SNR = L_signal - L_total_noise

# --- Final Output ---
print("This script calculates the Signal-to-Noise Ratio (SNR) at the people's location.")
print(f"1. Signal Level (L_S): {L_signal:.2f} dB")
print(f"2. Noise from train at the location (L_train): {L_train_at_people:.2f} dB")
print(f"3. Noise from construction at the location (L_const): {L_const_at_people:.2f} dB")
print(f"4. Total combined noise level (L_N): {L_total_noise:.2f} dB")
print("\nThe final equation for SNR is Signal Level - Total Noise Level:")
print(f"SNR = {L_signal:.2f} dB - {L_total_noise:.2f} dB = {SNR:.2f} dB")

<<<D>>>