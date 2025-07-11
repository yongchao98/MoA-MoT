import math

# Step 1: Define the signal level
L_signal = 75.0  # dB

# --- Calculate Total Noise Level ---

# Step 2: Calculate the train's noise level at the people's location
L_train_at_10m = 100.0  # dB
dist_train_to_people = 30.0  # meters
L_train_at_people = L_train_at_10m - 20 * math.log10(dist_train_to_people / 10.0)

# Step 3: Calculate the construction's noise level at the people's location
L_construction_at_20m = 115.0  # dB
dist_construction_to_people = 50.0  # meters
L_construction_at_people = L_construction_at_20m - 20 * math.log10(dist_construction_to_people / 20.0)

# Step 4: Combine the noise levels
# Convert dB to relative intensity
I_train_rel = 10**(L_train_at_people / 10.0)
I_construction_rel = 10**(L_construction_at_people / 10.0)

# Add intensities
I_noise_total_rel = I_train_rel + I_construction_rel

# Convert total intensity back to dB
L_noise_total = 10 * math.log10(I_noise_total_rel)

# Step 5: Calculate the final Signal-to-Noise Ratio (SNR)
snr = L_signal - L_noise_total

# --- Output the results ---
print("This script calculates the Signal-to-Noise Ratio (SNR) at the people's location.")
print(f"Signal Level (S): {L_signal:.2f} dB")
print(f"Calculated Train Noise at Location: {L_train_at_people:.2f} dB")
print(f"Calculated Construction Noise at Location: {L_construction_at_people:.2f} dB")
print(f"Total Combined Noise Level (N): {L_noise_total:.2f} dB")
print("\nThe final SNR is calculated as S - N:")
print(f"{L_signal:.2f} dB - {L_noise_total:.2f} dB = {snr:.2f} dB")
print("\nBased on the calculation, the closest answer is -32.")
