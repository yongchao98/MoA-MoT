import math

# Step 1: Define given values and the goal
L_signal = 75.0  # dB, Signal level at the location of the people

# Distances for the final SNR calculation
r_train_final = 30.0  # meters
r_const_final = 50.0  # meters

# Measurement 1 data
L_total_1 = 100.0  # dB
r_train_1 = 10.0   # meters
r_people_1 = 20.0  # meters

# Measurement 2 data
L_total_2 = 115.0  # dB
r_const_2 = 20.0   # meters
r_people_2 = 30.0  # meters

# Step 2: Assumption - The people's sound level at 1m is 75 dB
L_people_1m = 75.0

# Step 3: Calculate the people's contribution at the measurement points
# Sound level L at distance r is L(r) = L(r_ref) - 20 * log10(r / r_ref)
# Here r_ref = 1, so L(r) = L(1m) - 20 * log10(r)
L_people_at_M1 = L_people_1m - 20 * math.log10(r_people_1)
L_people_at_M2 = L_people_1m - 20 * math.log10(r_people_2)

# Step 4: Isolate the train and construction sound levels from the measurements
# Convert dB to normalized intensity (I/I_0 = 10^(L/10))
I_total_1_norm = 10**(L_total_1 / 10)
I_people_at_M1_norm = 10**(L_people_at_M1 / 10)
I_train_at_M1_norm = I_total_1_norm - I_people_at_M1_norm
L_train_at_M1 = 10 * math.log10(I_train_at_M1_norm) # dB at 10m

I_total_2_norm = 10**(L_total_2 / 10)
I_people_at_M2_norm = 10**(L_people_at_M2 / 10)
I_const_at_M2_norm = I_total_2_norm - I_people_at_M2_norm
L_const_at_M2 = 10 * math.log10(I_const_at_M2_norm) # dB at 20m

# Step 5: Calculate the noise levels at the final location (where the people are)
# L(r2) = L(r1) - 20 * log10(r2 / r1)
L_train_final = L_train_at_M1 - 20 * math.log10(r_train_final / r_train_1)
L_const_final = L_const_at_M2 - 20 * math.log10(r_const_final / r_const_2)

# Step 6: Combine the two noise sources to get the total noise level
I_train_final_norm = 10**(L_train_final / 10)
I_const_final_norm = 10**(L_const_final / 10)
I_noise_total_norm = I_train_final_norm + I_const_final_norm
L_noise = 10 * math.log10(I_noise_total_norm)

# Step 7: Calculate the final Signal-to-Noise Ratio (SNR)
SNR = L_signal - L_noise

# Print the final equation with all the numbers
print("Final Calculation:")
print(f"Signal Level (L_signal) = {L_signal:.2f} dB")
print(f"Noise Level (L_noise) = 10 * log10(10^({L_train_final:.2f}/10) + 10^({L_const_final:.2f}/10)) = {L_noise:.2f} dB")
print(f"Signal-to-Noise Ratio (SNR) = L_signal - L_noise")
print(f"SNR = {L_signal:.2f} dB - {L_noise:.2f} dB = {SNR:.2f} dB")
<<<D>>>