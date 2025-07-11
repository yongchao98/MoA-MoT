import math

# --- 1. Define Known and Assumed Parameters ---

# Number of lanes in one direction
N_L = 3
# Width of each lane in meters
W_L = 3.6
# Cross-slope in percent
S_percent = 1.75
# Pavement type described in the problem
pavement_type = "rough-textured asphalt"

# Manning's roughness coefficient 'n' for the given pavement type.
# This is a standard engineering value.
n = 0.016

# Design rainfall intensity 'i' in mm/hr.
# Since rainfall curves were not provided, a standard design value for
# hydroplaning analysis is assumed. 100 mm/hr is a common choice.
i_mm_hr = 100

# --- 2. Perform Calculations ---

# Calculate the total drainage path length 'L' in meters
L = N_L * W_L

# Convert the cross-slope 'S' from percent to a decimal
S = S_percent / 100.0

# Convert the rainfall intensity 'i' from mm/hr to m/s for the formula
# 1 hr = 3600 s, 1 m = 1000 mm
i_m_s = i_mm_hr / (3600 * 1000)

# --- 3. Apply the Kinematic Wave Formula ---

# Calculate the water film thickness 'd_m' in meters
# d_m = [ (n * L * i_m/s) / (S^0.5) ]^0.6
numerator = n * L * i_m_s
denominator = math.sqrt(S)
d_m = (numerator / denominator) ** 0.6

# Convert the final result 'd_mm' to millimeters
d_mm = d_m * 1000

# --- 4. Print the Results ---

print("--- Design Parameters ---")
print(f"Number of lanes: {N_L}")
print(f"Lane width: {W_L} m")
print(f"Total drainage path length (L): {L:.2f} m")
print(f"Cross-slope (S): {S_percent}% or {S}")
print(f"Pavement type: '{pavement_type}'")
print(f"Assumed Manning's n: {n}")
print(f"Assumed rainfall intensity (i): {i_mm_hr} mm/hr ({i_m_s:.3e} m/s)")
print("\n--- Calculation ---")
print("Using the Kinematic Wave equation: d_mm = 1000 * ((n * L * i_m/s) / (S**0.5))**0.6")
print("\nFinal equation with substituted values:")
# Using f-string formatting to display the equation with numbers
# The {:.3f} formats a float to 3 decimal places, {:.5f} to 5, etc.
# The :.3e formats to scientific notation with 3 decimal places.
equation_str = (
    f"d_mm = 1000 * (({n:.3f} * {L:.2f} * {i_m_s:.3e}) / ({S:.4f}**0.5))**0.6"
)
print(equation_str)

print(f"\nResulting water film thickness: {d_mm:.2f} mm")

# Return final answer in the requested format
# print(f"\n<<<{d_mm:.2f}>>>")