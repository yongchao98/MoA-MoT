import math

# Step 1: Define initial values from the problem.
L_signal = 75  # dB, Signal level at the location.
# Data point for Train: 100 dB at 10 meters.
Lt_measured = 100
r_train_measured = 10
# Data point for Construction: 115 dB at 20 meters.
Lc_measured = 115
r_const_measured = 20
# Distances of noise sources from the people's location.
r_train_final = 30
r_const_final = 50

# Step 2: Calculate the sound level from each noise source at the people's location.
# Using the formula: L2 = L1 - 20 * log10(r2 / r1).

# For the train:
Lt_final = Lt_measured - 20 * math.log10(r_train_final / r_train_measured)

# For the construction site:
Lc_final = Lc_measured - 20 * math.log10(r_const_final / r_const_measured)

# Step 3: Combine the two noise source levels to get the total noise level.
# Using the formula: L_total = 10 * log10(10^(L1/10) + 10^(L2/10)).
intensity_train = 10**(Lt_final / 10)
intensity_const = 10**(Lc_final / 10)
L_noise_total = 10 * math.log10(intensity_train + intensity_const)

# Step 4: Calculate the Signal-to-Noise Ratio (SNR).
SNR = L_signal - L_noise_total

# --- Output the results step-by-step ---
print(f"The signal level is given as {L_signal} dB.")
print("\n--- Noise Calculation ---")
print(f"1. The train's sound level at the location ({r_train_final}m) is calculated from its level at {r_train_measured}m:")
print(f"   Level = {Lt_measured} - 20 * log10({r_train_final}/{r_train_measured}) = {Lt_final:.2f} dB")

print(f"\n2. The construction's sound level at the location ({r_const_final}m) is calculated from its level at {r_const_measured}m:")
print(f"   Level = {Lc_measured} - 20 * log10({r_const_final}/{r_const_measured}) = {Lc_final:.2f} dB")

print("\n3. The total noise level is the combination of the train and construction noise:")
print(f"   Total Noise = 10 * log10(10^({Lt_final:.2f}/10) + 10^({Lc_final:.2f}/10)) = {L_noise_total:.2f} dB")

print("\n--- Final SNR Calculation ---")
print("The Signal-to-Noise Ratio (SNR) is the Signal Level minus the Total Noise Level.")
print(f"SNR = {L_signal} - {L_noise_total:.2f} = {SNR:.2f} dB")
<<<D>>>