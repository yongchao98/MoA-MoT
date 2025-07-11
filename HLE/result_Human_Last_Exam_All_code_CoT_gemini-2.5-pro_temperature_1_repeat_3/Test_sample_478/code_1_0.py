import math

# Step 1: Define the given signal level.
S_dB = 75.0

# Step 2: Define noise source levels and distances from measurements.
# We assume the people's contribution to these measurements is negligible.
L_train_ref = 100.0  # dB at 10m
d_train_ref = 10.0   # meters

L_const_ref = 115.0  # dB at 20m
d_const_ref = 20.0   # meters

# Step 3: Calculate the noise level from each source at the people's location.
d_train_final = 30.0 # meters
d_const_final = 50.0 # meters

# Sound level formula: L2 = L1 - 20 * log10(d2 / d1)
L_train_final = L_train_ref - 20 * math.log10(d_train_final / d_train_ref)
L_const_final = L_const_ref - 20 * math.log10(d_const_final / d_const_ref)

# Step 4: Convert the individual noise levels (dB) to relative intensities.
# Relative Intensity I_rel = 10^(L_dB / 10)
I_train_rel = 10**(L_train_final / 10)
I_const_rel = 10**(L_const_final / 10)

# Sum the intensities to find the total noise intensity.
I_total_noise_rel = I_train_rel + I_const_rel

# Step 5: Convert the total noise intensity back to a decibel level (N_dB).
N_dB = 10 * math.log10(I_total_noise_rel)

# Step 6: Calculate the Signal-to-Noise Ratio (SNR) in dB.
SNR_dB = S_dB - N_dB

# Step 7: Print the numbers in the final equation.
print(f"Final SNR Calculation:")
print(f"Signal Level: {S_dB:.2f} dB")
print(f"Total Noise Level: {N_dB:.2f} dB")
print(f"SNR = {S_dB:.2f} dB - {N_dB:.2f} dB = {SNR_dB:.2f} dB")