import math

# Step 1: Define the given parameters and assumptions
lanes_per_direction = 3
lane_width_m = 3.6  # in meters
cross_slope_percent = 1.75

# Manning's roughness coefficient for rough-textured asphalt (standard assumption)
n_manning = 0.015

# Design storm parameters for hydroplaning analysis (standard assumption)
# A 10-year, 5-minute storm is common for this type of critical check.
storm_duration_min = 5

# Generic IDF (Intensity-Duration-Frequency) formula: i = a / (t + b)
# Using common constants for i in mm/hr and t in minutes.
idf_a = 3048
idf_b = 20

# Step 2: Calculate intermediate values
# Calculate the total flow path length (L) in meters
flow_path_length_m = lanes_per_direction * lane_width_m

# Convert the cross-slope from percent to decimal
cross_slope_decimal = cross_slope_percent / 100

# Calculate the design rainfall intensity (i) in mm/hr
rainfall_intensity_mmhr = idf_a / (storm_duration_min + idf_b)

# Step 3: Calculate the water film thickness using the FHWA formula
# The formula for Tw in mm is: Tw = 1.22 * (n * L)^0.6 * i^0.4 * S^(-0.2)
# where L is in meters, i is in mm/hr, and S is a decimal.

# Calculate each part of the formula
term1 = n_manning * flow_path_length_m
term2 = rainfall_intensity_mmhr
term3 = cross_slope_decimal

# Calculate the final water film thickness in mm
water_film_thickness_mm = 1.22 * math.pow(term1, 0.6) * math.pow(term2, 0.4) * math.pow(term3, -0.2)

# Step 4: Print the results and the final equation
print("--- Design Water Film Thickness Calculation ---")
print(f"Flow Path Length (L): {lanes_per_direction} lanes * {lane_width_m} m/lane = {flow_path_length_m:.1f} m")
print(f"Pavement Cross-Slope (S): {cross_slope_percent}% = {cross_slope_decimal}")
print(f"Manning's Roughness (n): {n_manning} (Assumed for rough-textured asphalt)")
print(f"Design Rainfall Intensity (i): {rainfall_intensity_mmhr:.2f} mm/hr (For a 10-year, 5-minute storm)")
print("\n--- Final Equation ---")
print("The water film thickness (T_w) is calculated using the formula:")
print("T_w(mm) = 1.22 * (n * L)^0.6 * i^0.4 * S^(-0.2)\n")

# Print the equation with the calculated values
print("Plugging in the values:")
print(f"T_w(mm) = 1.22 * ({n_manning} * {flow_path_length_m:.1f})^0.6 * {rainfall_intensity_mmhr:.2f}^0.4 * {cross_slope_decimal}^(-0.2)")

# Print the final calculated result
print("\n--- Result ---")
print(f"The design water film thickness is: {water_film_thickness_mm:.2f} mm")
print("<<<6.80>>>")