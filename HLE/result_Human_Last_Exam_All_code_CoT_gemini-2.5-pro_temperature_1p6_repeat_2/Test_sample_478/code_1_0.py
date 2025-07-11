import math

def add_spl(levels):
    """Combines multiple sound pressure levels (in dB)."""
    total_intensity_ratio = sum(10**(L / 10) for L in levels)
    if total_intensity_ratio <= 0:
        return -float('inf')
    return 10 * math.log10(total_intensity_ratio)

def subtract_spl(total_level, known_level):
    """Subtracts a known sound source from a total SPL (in dB)."""
    if known_level >= total_level:
        # If the 'known' sound is louder than the total, it implies the other source is negligible.
        # This can happen due to small floating point inaccuracies. In reality the other source SPL is very low.
        # We handle this by calculating based on the numbers provided.
        # But if total is less, it's a very tiny difference that is being subtracted.
        # If we return -inf, it would cause issues. Instead lets assume the resulting level is very small but not infinitely small
        # A more robust check:
        total_intensity = 10**(total_level / 10)
        known_intensity = 10**(known_level / 10)
        if known_intensity >= total_intensity:
           return -float('inf') # Should mathematically not be the primary source, practically negligible
        
    resulting_intensity_ratio = 10**(total_level / 10) - 10**(known_level / 10)
    if resulting_intensity_ratio <= 0:
        return -float('inf') 
    return 10 * math.log10(resulting_intensity_ratio)

def project_spl(L1, r1, r2):
    """Calculates SPL at distance r2, given SPL L1 at distance r1."""
    if r1 <= 0 or r2 <= 0:
        raise ValueError("Distances must be positive")
    return L1 - 20 * math.log10(r2 / r1)

# --- Problem Data ---
# Signal Level at the point of interest
L_signal = 75.0  # dB

# Measurement 1: Train + People
d1_train = 10.0   # meters
d1_people = 20.0  # meters
L_total1 = 100.0  # dB

# Measurement 2: Construction + People
d2_constr = 20.0  # meters
d2_people = 30.0  # meters
L_total2 = 115.0  # dB

# Final Location (where people are, i.e., the receiver's location)
d_final_train = 30.0   # meters
d_final_constr = 50.0  # meters

# --- Calculations ---
# Assumption: The 75 dB of the people is their SPL at a reference distance of 1m.
L_people_1m = 75.0

# Step 1: Analyze Measurement 1 to find Train's source strength
# SPL of people at measurement point 1 (20m away)
L_people_at_M1 = project_spl(L_people_1m, 1, d1_people)
# SPL of train alone at measurement point 1 (10m away)
L_train_at_M1 = subtract_spl(L_total1, L_people_at_M1)
# Source strength of train (projected back to 1m)
L_train_1m = project_spl(L_train_at_M1, d1_train, 1)

# Step 2: Analyze Measurement 2 to find Construction's source strength
# SPL of people at measurement point 2 (30m away)
L_people_at_M2 = project_spl(L_people_1m, 1, d2_people)
# SPL of construction alone at measurement point 2 (20m away)
L_constr_at_M2 = subtract_spl(L_total2, L_people_at_M2)
# Source strength of construction (projected back to 1m)
L_constr_1m = project_spl(L_constr_at_M2, d2_constr, 1)

# Step 3: Calculate noise SPL from each source at the final location
# SPL of train at the people's location (30m away)
L_train_final = project_spl(L_train_1m, 1, d_final_train)
# SPL of construction at the people's location (50m away)
L_constr_final = project_spl(L_constr_1m, 1, d_final_constr)

# Step 4: Calculate total noise level at the final location
L_noise_total = add_spl([L_train_final, L_constr_final])

# Step 5: Calculate Signal-to-Noise Ratio (SNR)
snr = L_signal - L_noise_total

# --- Output ---
print("Calculating Signal-to-Noise Ratio (SNR)\n")
print("--- Final Equation ---")
print("SNR = Signal Level - Total Noise Level\n")
print(f"Signal Level = {L_signal:.2f} dB")
print(f"Noise from Train at location = {L_train_final:.2f} dB")
print(f"Noise from Construction at location = {L_constr_final:.2f} dB")
print(f"Total Noise Level (Train + Construction) = {L_noise_total:.2f} dB\n")
print("--- Result ---")
print(f"SNR = {L_signal:.2f} - {L_noise_total:.2f}")
print(f"SNR = {snr:.2f} dB")
