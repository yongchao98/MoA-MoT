import math

# --- Step 1: Define given and assumed parameters ---
# Given parameters from the problem description
num_lanes = 3
lane_width_m = 3.6  # meters
cross_slope_percent = 1.75  # %

# Assumed parameters based on standard engineering practice
# Manning's n for rough-textured asphalt pavement
manning_n = 0.016
# Design rainfall intensity for hydroplaning analysis
rainfall_intensity_mm_hr = 100  # mm/hr

# --- Step 2: Calculate derived parameters for the formula ---
# Drainage path length (L) is the total width of the three lanes
L_m = num_lanes * lane_width_m

# Convert cross-slope from percent to a decimal value (m/m)
S_decimal = cross_slope_percent / 100.0

# Convert rainfall intensity from mm/hr to m/s for the SI unit formula
# 1 mm/hr = 1 / (1000 mm/m) / (3600 s/hr) = 1 / 3,600,000 m/s
i_m_s = rainfall_intensity_mm_hr / 3600000.0

# --- Step 3: Calculate the water film thickness ---
# Using the Kinematic Wave formula for overland sheet flow
# T_w = [ (n * i * L) / (S^0.5) ]^0.6
numerator = manning_n * i_m_s * L_m
denominator = math.sqrt(S_decimal)
Tw_m = (numerator / denominator) ** 0.6

# Convert final result from meters to millimeters
Tw_mm = Tw_m * 1000

# --- Step 4: Print the final result and equation ---
print("This script calculates the design water film thickness on a pavement surface.")
print("-" * 70)
print(f"Given and Assumed Parameters:")
print(f"  - Drainage Path Length (L): {L_m:.2f} m ({num_lanes} lanes @ {lane_width_m} m)")
print(f"  - Pavement Cross-Slope (S): {S_decimal:.4f} m/m ({cross_slope_percent}%)")
print(f"  - Manning's Roughness (n): {manning_n} (for rough-textured asphalt)")
print(f"  - Design Rainfall Intensity (i): {rainfall_intensity_mm_hr} mm/hr ({i_m_s:.6f} m/s)")
print("-" * 70)

print("The final calculation using the Kinematic Wave equation is:")
# Print the equation with the variables plugged in, as requested.
# The calculation for Tw_m is performed, and the final result is multiplied by 1000 for mm.
print(f"T_w (mm) = 1000 * [ (n * i_m_s * L_m) / (S_decimal**0.5) ]**0.6")
print(f"T_w (mm) = 1000 * [ ({manning_n} * {i_m_s:.6f} * {L_m:.2f}) / ({S_decimal:.4f}**0.5) ]**0.6")
print(f"T_w (mm) = {Tw_mm:.2f}")
print("-" * 70)
print(f"The design water film thickness at the outer edge is {Tw_mm:.2f} mm.")
