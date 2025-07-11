import math

# Step 1: Define the signal level
S = 75  # dB, Signal level at the people's location

# Step 2 & 3: Define noise source levels at their measurement points
# Noise from the train
L1_train = 100  # dB
r1_train = 10   # meters

# Noise from the construction
L1_constr = 115  # dB
r1_constr = 20   # meters

# Define the distances of the noise sources from the people's location
r2_train = 30   # meters
r2_constr = 50  # meters

# Step 4: Calculate the noise level from each source at the people's location
# Using the formula: L2 = L1 - 20 * log10(r2 / r1)

# Noise from train at the people's location
L2_train = L1_train - 20 * math.log10(r2_train / r1_train)

# Noise from construction at the people's location
L2_constr = L1_constr - 20 * math.log10(r2_constr / r1_constr)

print(f"Noise level from train at the location = {L2_train:.2f} dB")
print(f"Noise level from construction at the location = {L2_constr:.2f} dB")

# Step 5: Combine the two noise levels
# Using the formula: L_total = 10 * log10(10^(L1/10) + 10^(L2/10))
intensity_train = 10**(L2_train / 10)
intensity_constr = 10**(L2_constr / 10)
N_total = 10 * math.log10(intensity_train + intensity_constr)

print(f"Total noise level at the location = {N_total:.2f} dB")

# Step 6: Calculate the Signal-to-Noise Ratio (SNR)
SNR = S - N_total

# Output the final equation and result
print("\nFinal Calculation:")
print(f"SNR (dB) = Signal (dB) - Total Noise (dB)")
print(f"SNR = {S} - {N_total:.2f} = {SNR:.2f} dB")

# Step 7: Compare with options and find the closest one.
# The calculated SNR is approximately -32.14 dB.
# Answer Choices: A. -36.6, B. -26.13, C. -35.41, D. -32, E. -15
# The closest answer is D.
