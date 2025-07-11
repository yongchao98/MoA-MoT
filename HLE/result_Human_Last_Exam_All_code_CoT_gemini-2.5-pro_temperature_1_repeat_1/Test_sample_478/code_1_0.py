import math

# Step 1: Define the known values from the problem statement.
signal_level = 75.0  # dB, at the people's location

# Data for the train
train_level_ref = 100.0  # dB
train_dist_ref = 10.0  # meters
train_dist_at_people = 30.0  # meters

# Data for the construction
construction_level_ref = 115.0  # dB
construction_dist_ref = 20.0  # meters
construction_dist_at_people = 50.0  # meters

# Step 2: Calculate the noise level from each source at the people's location.

# Noise from the train at the people's location
train_level_at_people = train_level_ref - 20 * math.log10(train_dist_at_people / train_dist_ref)

# Noise from the construction at the people's location
construction_level_at_people = construction_level_ref - 20 * math.log10(construction_dist_at_people / construction_dist_ref)

# Step 3: Combine the noise levels.
# Convert dB to intensity-like values, sum them, and convert back to dB.
intensity_train = 10**(train_level_at_people / 10)
intensity_construction = 10**(construction_level_at_people / 10)
total_noise_intensity = intensity_train + intensity_construction
total_noise_level = 10 * math.log10(total_noise_intensity)

# Step 4: Calculate the Signal-to-Noise Ratio (SNR).
snr = signal_level - total_noise_level

# Print the results step-by-step
print(f"Signal Level: {signal_level:.2f} dB")
print(f"Train noise level at the people's location: {train_level_at_people:.2f} dB")
print(f"Construction noise level at the people's location: {construction_level_at_people:.2f} dB")
print(f"Total combined noise level: {total_noise_level:.2f} dB")
print("\nFinal Calculation:")
print(f"Signal-to-Noise Ratio (SNR) = Signal Level - Total Noise Level")
print(f"SNR = {signal_level:.2f} dB - {total_noise_level:.2f} dB = {snr:.2f} dB")

<<<D>>>