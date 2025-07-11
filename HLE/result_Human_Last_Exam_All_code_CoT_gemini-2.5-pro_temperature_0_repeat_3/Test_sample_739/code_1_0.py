import math

# --- Step 1: Define given parameters and assumptions ---

# Given parameters from the problem
num_lanes = 3
lane_width_m = 3.6  # meters
cross_slope_percent = 1.75 # %

# Assumptions for missing data
# Manning's roughness coefficient 'n' for rough-textured asphalt (sheet flow)
n = 0.016
# Design rainfall intensity 'i' in mm/hr. A value typical for hydroplaning analysis.
i_mm_hr = 125.0

# --- Step 2: Prepare variables for the formula ---

# Calculate total flow path length in meters
L_m = num_lanes * lane_width_m

# Convert cross-slope from percent to decimal
S = cross_slope_percent / 100.0

# The HEC-22 formula is in US Customary Units, so we convert our inputs.
# Conversion factors
M_TO_FT = 3.28084
MM_TO_IN = 1.0 / 25.4

# Convert length to feet
L_ft = L_m * M_TO_FT

# Convert rainfall intensity to inches per hour
i_in_hr = i_mm_hr * MM_TO_IN

# --- Step 3: Calculate water film thickness using the HEC-22 formula ---

# HEC-22 Formula: d_in = 0.007 * (n * L_ft)^0.8 * i_in_hr^0.6 / S^0.4
# Calculate the water depth in inches
d_in = 0.007 * math.pow(n * L_ft, 0.8) * math.pow(i_in_hr, 0.6) / math.pow(S, 0.4)

# Convert the final result from inches to millimeters
d_mm = d_in * 25.4

# --- Step 4: Print the results ---

print("Calculation of Water Film Thickness (d_w)")
print("-" * 40)
print(f"Inputs and Assumptions:")
print(f"  - Flow Path Length (L): {num_lanes} lanes * {lane_width_m} m/lane = {L_m:.2f} m ({L_ft:.2f} ft)")
print(f"  - Cross-Slope (S): {cross_slope_percent}% = {S}")
print(f"  - Assumed Manning's n: {n}")
print(f"  - Assumed Rainfall Intensity (i): {i_mm_hr} mm/hr ({i_in_hr:.2f} in/hr)")
print("-" * 40)

# Print the final equation with all the numbers plugged in
# The formula calculates depth in inches, which is then converted to mm
print("Formula (HEC-22): d_w = (0.007 * (n * L_ft)^0.8 * i_in_hr^0.6 / S^0.4) * 25.4")
print("\nFinal Equation:")
print(f"d_w = (0.007 * ({n} * {L_ft:.2f})^0.8 * {i_in_hr:.2f}^0.6 / {S}^0.4) * 25.4")
print(f"\nResult:")
print(f"The design water film thickness is {d_mm:.2f} mm.")
<<<1.57>>>