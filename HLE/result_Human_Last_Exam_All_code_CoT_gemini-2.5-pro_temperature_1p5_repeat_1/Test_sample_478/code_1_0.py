import math

def add_db(levels):
    """Correctly adds a list of sound levels in decibels."""
    if not levels:
        return -float('inf')
    # Sum the corresponding intensities
    total_intensity = sum(10**(level / 10.0) for level in levels)
    # Convert the total intensity back to decibels
    return 10 * math.log10(total_intensity)

def subtract_db(total_level, known_level):
    """
    Subtracts a known sound level from a total sound level.
    This is used to find the level of an unknown source when the combined
    level and the level of another source are known.
    """
    # If the known sound is louder than or equal to the total, the math is invalid.
    if total_level <= known_level:
        # The contribution of the unknown source is negligible or zero.
        return -float('inf') # Represents zero intensity
    
    # Convert levels to intensities
    total_intensity = 10**(total_level / 10.0)
    known_intensity = 10**(known_level / 10.0)
    
    # Subtract intensities
    result_intensity = total_intensity - known_intensity
    
    # Convert the resulting intensity back to decibels
    return 10 * math.log10(result_intensity)

def project_db(level_at_r1, r1, r2):
    """
    Calculates the sound level at distance r2, given the sound level at distance r1,
    assuming a point source (inverse square law).
    """
    if r1 <= 0 or r2 <= 0:
        raise ValueError("Distances must be positive.")
    # The change in level is 20 * log10(r1 / r2)
    return level_at_r1 + 20 * math.log10(r1 / r2)

# --- Problem Data ---
# Signal level at the people's location
L_signal = 75.0

# Assumption: The signal level is based on a source strength of 75 dB at 1 meter.
L_p_1m = 75.0

# Measurement 1 (Train + People)
L_M1_total = 100.0
d_t_M1 = 10.0
d_p_M1 = 20.0

# Measurement 2 (Construction + People)
L_M2_total = 115.0
d_c_M2 = 20.0
d_p_M2 = 30.0

# Distances from noise sources to the people's location
d_t_p = 30.0
d_c_p = 50.0

# --- Calculations ---

# 1. Determine the sound level of the TRAIN at the people's location
# 1a. Find the people's sound level at the first measurement point (M1).
L_p_at_M1 = project_db(L_p_1m, 1.0, d_p_M1)
# 1b. Find the train's sound level at M1 by removing the people's contribution.
L_t_at_M1 = subtract_db(L_M1_total, L_p_at_M1)
# 1c. Project the train's sound from M1 to the people's location.
L_t_at_p = project_db(L_t_at_M1, d_t_M1, d_t_p)

# 2. Determine the sound level of the CONSTRUCTION at the people's location
# 2a. Find the people's sound level at the second measurement point (M2).
L_p_at_M2 = project_db(L_p_1m, 1.0, d_p_M2)
# 2b. Find the construction's sound level at M2 by removing the people's contribution.
L_c_at_M2 = subtract_db(L_M2_total, L_p_at_M2)
# 2c. Project the construction's sound from M2 to the people's location.
L_c_at_p = project_db(L_c_at_M2, d_c_M2, d_c_p)

# 3. Calculate the TOTAL NOISE level at the people's location by adding the two noise sources.
L_noise_total = add_db([L_t_at_p, L_c_at_p])

# 4. Calculate the Signal-to-Noise Ratio (SNR).
snr = L_signal - L_noise_total

# --- Output the Final Answer and Its Components ---
print("--- Calculating Signal-to-Noise Ratio (SNR) ---")
print(f"The defined signal level is: {L_signal:.2f} dB")
print(f"Noise from Train at the location is calculated to be: {L_t_at_p:.2f} dB")
print(f"Noise from Construction at the location is calculated to be: {L_c_at_p:.2f} dB")
print(f"The total combined noise level is: {L_noise_total:.2f} dB")
print("\nThe final equation is: SNR = Signal Level - Total Noise Level")
print(f"SNR = {L_signal:.2f} dB - {L_noise_total:.2f} dB")
print(f"The final Signal-to-Noise Ratio (SNR) is: {snr:.2f} dB")
