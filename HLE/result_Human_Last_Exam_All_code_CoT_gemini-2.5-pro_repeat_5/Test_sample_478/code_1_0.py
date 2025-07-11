import math

# --- Given Information ---
# Signal
L_signal = 75  # dB, sound level of people at their location

# Noise source distances from the people
d_train_to_people = 30  # meters
d_construction_to_people = 50  # meters

# Measurement 1 (Train + People)
L_combo1 = 100  # dB
d_train_to_m1 = 10  # meters
d_people_to_m1 = 20 # meters

# Measurement 2 (Construction + People)
L_combo2 = 115  # dB
d_construction_to_m2 = 20  # meters
d_people_to_m2 = 30 # meters

# --- Step 1: Analyze the Train Noise ---
# The contribution of the people's noise (around 49 dB) to the 100 dB measurement is negligible.
# So, we can approximate the train's sound level at measurement point 1.
L_train_at_m1 = L_combo1

# Calculate the train's sound level at the people's location (30m)
# using the formula L2 = L1 - 20*log10(r2/r1)
L_train_at_people = L_train_at_m1 - 20 * math.log10(d_train_to_people / d_train_to_m1)

# --- Step 2: Analyze the Construction Noise ---
# The contribution of the people's noise (around 45 dB) to the 115 dB measurement is negligible.
# So, we can approximate the construction's sound level at measurement point 2.
L_construction_at_m2 = L_combo2

# Calculate the construction's sound level at the people's location (50m)
L_construction_at_people = L_construction_at_m2 - 20 * math.log10(d_construction_to_people / d_construction_to_m2)

# --- Step 3: Calculate Total Noise at the People's Location ---
# Convert dB levels of noise sources to relative intensities
I_train_rel = 10**(L_train_at_people / 10)
I_construction_rel = 10**(L_construction_at_people / 10)

# Add the intensities to get the total noise intensity
I_noise_total_rel = I_train_rel + I_construction_rel

# Convert the total intensity back to dB
L_noise_total = 10 * math.log10(I_noise_total_rel)

# --- Step 4: Calculate the Signal-to-Noise Ratio (SNR) ---
snr = L_signal - L_noise_total

# --- Final Output ---
print("--- Calculating the Signal-to-Noise Ratio (SNR) ---")
print(f"Signal Level (people talking): {L_signal} dB")
print("\nCalculating Noise Components at the People's Location:")
print(f"1. Train Noise Level = {L_train_at_people:.2f} dB")
print(f"2. Construction Noise Level = {L_construction_at_people:.2f} dB")
print("\nCombining noise sources to find total noise level:")
# Using f-string formatting to display the equation with calculated values
equation_noise = f"L_noise = 10 * log10(10^({L_train_at_people:.2f}/10) + 10^({L_construction_at_people:.2f}/10))"
print(f"{equation_noise} = {L_noise_total:.2f} dB")

print("\nCalculating the final SNR:")
equation_snr = f"SNR = Signal Level - Total Noise Level"
print(equation_snr)
final_calculation = f"SNR = {L_signal} dB - {L_noise_total:.2f} dB = {snr:.2f} dB"
print(final_calculation)
print("\nThe final Signal-to-Noise Ratio is approximately -32 dB.")
<<<D>>>