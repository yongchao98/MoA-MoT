import math

# Step 1: Define the given and assumed parameters
# Road Geometry
num_lanes_per_direction = 3
lane_width_m = 3.6  # in meters
cross_slope_percent = 1.75

# Hydraulic Parameters (Assumed based on standard engineering practice)
# Manning's n for rough-textured asphalt pavement
manning_n = 0.016
# Design rainfall intensity (i) in mm/hr for hydroplaning analysis
rainfall_intensity_mmhr = 150.0

# Step 2: Calculate intermediate values
# L: Total drainage path length in meters
drainage_path_length_L = num_lanes_per_direction * lane_width_m
# S: Cross-slope as a decimal (m/m)
cross_slope_S = cross_slope_percent / 100.0

# Step 3: Calculate the water film thickness using the kinematic wave formula
# Formula: Tw = C * (n * L)^0.6 * i^0.6 * S^-0.3
# Where C is a constant (0.3051) to yield Tw in mm from L in m and i in mm/hr.
conversion_constant = 0.3051

# Perform the calculation
water_film_thickness_mm = (
    conversion_constant *
    math.pow(manning_n * drainage_path_length_L, 0.6) *
    math.pow(rainfall_intensity_mmhr, 0.6) *
    math.pow(cross_slope_S, -0.3)
)

# Step 4: Print the results, showing the full equation with values
print("--- Design Water Film Thickness Calculation ---")
print("\nParameters:")
print(f"  - Drainage Path Length (L): {num_lanes_per_direction} lanes * {lane_width_m} m/lane = {drainage_path_length_L:.1f} m")
print(f"  - Cross-Slope (S): {cross_slope_percent}% = {cross_slope_S}")
print(f"  - Manning's n (for rough asphalt): {manning_n} (Assumed)")
print(f"  - Design Rainfall Intensity (i): {rainfall_intensity_mmhr:.1f} mm/hr (Assumed)")

print("\nFormula:")
print("  Tw = C * (n * L)^0.6 * i^0.6 * S^-0.3")

print("\nCalculation:")
print(f"  Tw = {conversion_constant:.4f} * ({manning_n} * {drainage_path_length_L:.1f})^0.6 * {rainfall_intensity_mmhr:.1f}^0.6 * {cross_slope_S}^-0.3")

# Calculating each part of the equation for clarity in the output
part1 = math.pow(manning_n * drainage_path_length_L, 0.6)
part2 = math.pow(rainfall_intensity_mmhr, 0.6)
part3 = math.pow(cross_slope_S, -0.3)
print(f"  Tw = {conversion_constant:.4f} * ({part1:.4f}) * ({part2:.4f}) * ({part3:.4f})")

print("\nResult:")
print(f"The design water film thickness is {water_film_thickness_mm:.2f} mm.")
<<<6.01>>>