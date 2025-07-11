import math

# Step 1: Define the known values from the problem statement.
S_dB = 75.0  # Signal level in dB

# Train noise parameters
L_t_10m = 100.0  # Train noise level at 10m
r_t_ref = 10.0   # Reference distance for train noise
r_t_final = 30.0 # Final distance of train from people

# Construction noise parameters
L_c_20m = 115.0  # Construction noise level at 20m
r_c_ref = 20.0   # Reference distance for construction noise
r_c_final = 50.0 # Final distance of construction from people

# Step 2: Calculate the noise level from the train at the people's location.
# Formula: L2 = L1 + 20 * log10(r1 / r2)
L_t_final = L_t_10m + 20 * math.log10(r_t_ref / r_t_final)

# Step 3: Calculate the noise level from the construction at the people's location.
L_c_final = L_c_20m + 20 * math.log10(r_c_ref / r_c_final)

# Step 4: Combine the noise sources.
# Convert dB to intensity (proportional value), add them, then convert back to dB.
# I = 10^(L/10)
I_t_final = 10**(L_t_final / 10)
I_c_final = 10**(L_c_final / 10)

# Total noise intensity is the sum of individual intensities
I_total_noise = I_t_final + I_c_final

# Convert total intensity back to dB for the total noise level N
# N = 10 * log10(I_total)
N_dB = 10 * math.log10(I_total_noise)

# Step 5: Calculate the Signal-to-Noise Ratio (SNR).
# SNR = Signal (dB) - Noise (dB)
SNR = S_dB - N_dB

# --- Output the results ---
print("--- Calculation Steps ---")
print(f"Signal Level (S): {S_dB} dB")
print("\nCalculating Total Noise Level (N):")
print(f"1. Noise from train at people's location ({r_t_final}m) = {L_t_10m} + 20*log10({r_t_ref}/{r_t_final}) = {L_t_final:.2f} dB")
print(f"2. Noise from construction at people's location ({r_c_final}m) = {L_c_20m} + 20*log10({r_c_ref}/{r_c_final}) = {L_c_final:.2f} dB")
print(f"3. Total noise level (N) = 10*log10( 10^({L_t_final:.2f}/10) + 10^({L_c_final:.2f}/10) ) = {N_dB:.2f} dB")
print("\n--- Final Calculation ---")
print(f"Signal-to-Noise Ratio (SNR) = Signal Level (S) - Total Noise Level (N)")
print(f"SNR = {S_dB} dB - {N_dB:.2f} dB = {SNR:.2f} dB")
