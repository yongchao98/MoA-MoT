import math

# Step 1: Define given and assumed parameters
num_lanes = 3
lane_width_m = 3.6  # in meters
cross_slope_percent = 1.75  # in percent

# Assumptions for missing data
# A standard Manning's n for rough-textured asphalt pavement
manning_n = 0.016
# A design rainfall intensity for a 10-year, 5-minute storm (typical for arterial roads)
rainfall_intensity_mmhr = 150.0 # in mm/hr

print("--- Design Parameters ---")
print(f"Number of lanes: {num_lanes}")
print(f"Lane width: {lane_width_m} m")
print(f"Cross-slope: {cross_slope_percent}%")
print(f"Assumed Manning's n for rough-textured asphalt: {manning_n}")
print(f"Assumed design rainfall intensity (I): {rainfall_intensity_mmhr} mm/hr\n")


# Step 2: Calculate derived parameters for the formula
# Drainage path length (L) in meters
L = num_lanes * lane_width_m
# Cross-slope (S) as a decimal
S = cross_slope_percent / 100

# Step 3: Apply the metric FHWA formula for water film thickness (Tw)
# Tw (mm) = 0.0847 * (n * L)^0.8 * I^0.6 * S^-0.4
term1 = math.pow(manning_n * L, 0.8)
term2 = math.pow(rainfall_intensity_mmhr, 0.6)
term3 = math.pow(S, -0.4)
Tw_mm = 0.0847 * term1 * term2 * term3

# Step 4: Print the final calculation and result
print("--- Calculation ---")
print("Formula: Tw (mm) = 0.0847 * (n * L)^0.8 * I^0.6 * S^-0.4")
print("Substituting values:")
# The problem asks to output each number in the final equation.
# To make it clear, we show the combined value of (n*L)
n_L_val = manning_n * L
print(f"Tw = 0.0847 * ({manning_n} * {L})^0.8 * {rainfall_intensity_mmhr}^0.6 * {S}^-0.4")
print(f"Tw = 0.0847 * ({n_L_val:.4f})^0.8 * {rainfall_intensity_mmhr}^0.6 * {S}^-0.4")
print(f"Tw = 0.0847 * {term1:.4f} * {term2:.4f} * {term3:.4f}")
print(f"\nFinal Result:")
print(f"The design water film thickness is {Tw_mm:.2f} mm.")
