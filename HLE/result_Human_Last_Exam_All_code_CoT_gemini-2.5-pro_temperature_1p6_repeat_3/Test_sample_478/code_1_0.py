import math

# 1. Define Signal Level and reference noise measurements.
# Signal level (S) is the sound of the people talking.
signal_db = 75.0

# Reference level for the train.
level_train_ref = 100.0
dist_train_ref = 10.0

# Reference level for the construction site.
level_const_ref = 115.0
dist_const_ref = 20.0

# 2. Define distances from noise sources to the people.
dist_train_to_people = 30.0
dist_const_to_people = 50.0

# 3. Calculate the noise level from the train at the people's location.
# Using the formula: L2 = L1 - 20 * log10(r2 / r1)
level_train_at_people = level_train_ref - 20 * math.log10(dist_train_to_people / dist_train_ref)

# 4. Calculate the noise level from the construction site at the people's location.
level_const_at_people = level_const_ref - 20 * math.log10(dist_const_to_people / dist_const_ref)

# 5. Combine the two noise levels to get the total noise.
# Convert dB to linear intensity, sum them up, and convert back to dB.
intensity_train = 10**(level_train_at_people / 10)
intensity_const = 10**(level_const_at_people / 10)
total_intensity_noise = intensity_train + intensity_const
total_noise_db = 10 * math.log10(total_intensity_noise)

# 6. Calculate the final Signal-to-Noise Ratio (SNR).
snr = signal_db - total_noise_db

# Print the results, including the final equation as requested.
print(f"Signal Level (S): {signal_db:.2f} dB")
print(f"Train Noise at Location (N1): {level_train_at_people:.2f} dB")
print(f"Construction Noise at Location (N2): {level_const_at_people:.2f} dB")
print(f"Total Noise Level (N_total): {total_noise_db:.2f} dB")
print("\nFinal SNR Calculation:")
print(f"SNR = Signal (dB) - Total Noise (dB)")
print(f"{signal_db:.2f} - {total_noise_db:.2f} = {snr:.2f} dB")
