import math

# Step 1: Define the signal level
L_signal = 75  # dB

# --- Calculate Noise ---

# Step 2: Calculate the noise level from the train at the people's location
L_train_at_10m = 100
distance_train_to_people = 30
L_train_at_people = L_train_at_10m - 20 * math.log10(distance_train_to_people / 10)

# Step 3: Calculate the noise level from the construction at the people's location
L_construction_at_20m = 115
distance_construction_to_people = 50
L_construction_at_people = L_construction_at_20m - 20 * math.log10(distance_construction_to_people / 20)

# Step 4: Combine the two noise sources to get total noise level
# Convert dB to relative intensity, sum them, and convert back to dB
I_train_relative = 10**(L_train_at_people / 10)
I_construction_relative = 10**(L_construction_at_people / 10)
L_noise_total = 10 * math.log10(I_train_relative + I_construction_relative)

# Step 5: Calculate the final Signal-to-Noise Ratio (SNR)
SNR = L_signal - L_noise_total

# --- Print the results ---
print(f"Signal Level: {L_signal:.2f} dB")
print(f"Noise from train at the location: {L_train_at_people:.2f} dB")
print(f"Noise from construction at the location: {L_construction_at_people:.2f} dB")
print(f"Total combined noise level: {L_noise_total:.2f} dB")
print("-" * 30)
print("The final equation for SNR is: Signal Level - Total Noise Level")
print(f"SNR (dB) = {L_signal} - {L_noise_total:.2f} = {SNR:.2f}")