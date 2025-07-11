import math

# Step 1: Define the known values from the problem statement.
L_signal = 75  # dB, the sound level of the people's conversation.

# Train data
L_train_measured = 100      # dB, the sound level measured
r_train_measured = 10       # meters, distance for the measurement
r_train_to_people = 30      # meters, distance from the train to the people

# Construction data
L_const_measured = 115      # dB, the sound level measured
r_const_measured = 20       # meters, distance for the measurement
r_const_to_people = 50      # meters, distance from the construction to the people

# Step 2: Calculate the train's noise level at the people's location.
# The formula for sound level change with distance is L2 = L1 + 20 * log10(r1 / r2)
L_train_at_people = L_train_measured + 20 * math.log10(r_train_measured / r_train_to_people)
print(f"Calculated noise level from the train at the people's location: {L_train_at_people:.2f} dB")

# Step 3: Calculate the construction's noise level at the people's location.
L_const_at_people = L_const_measured + 20 * math.log10(r_const_measured / r_const_to_people)
print(f"Calculated noise level from construction at the people's location: {L_const_at_people:.2f} dB")

# Step 4: Combine the two noise sources to find the total noise level.
# We must convert from dB to intensity, sum the intensities, then convert back to dB.
# Formula: L_total = 10 * log10(10^(L1/10) + 10^(L2/10))
intensity_train = 10**(L_train_at_people / 10)
intensity_const = 10**(L_const_at_people / 10)
L_total_noise = 10 * math.log10(intensity_train + intensity_const)
print(f"Total combined noise level at the people's location: {L_total_noise:.2f} dB")

# Step 5: Calculate the Signal-to-Noise Ratio (SNR) in dB.
# SNR (dB) = L_signal - L_noise
SNR = L_signal - L_total_noise

# Print the final equation and the result.
print("\n--- Final Answer Calculation ---")
print("The Signal-to-Noise Ratio (SNR) is calculated by subtracting the total noise level from the signal level.")
print(f"SNR (dB) = {L_signal} dB - {L_total_noise:.2f} dB")
print(f"The final Signal-to-Noise ratio is: {SNR:.2f} dB")
<<<D>>>