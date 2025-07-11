import math

# Step 1: Define the given and assumed parameters
# --- Geometric Parameters ---
num_lanes = 3
lane_width_m = 3.6  # in meters
cross_slope_percent = 1.75

# --- Hydraulic Parameters (Assumptions) ---
# Manning's roughness coefficient 'n' for rough-textured asphalt pavement.
# This is a standard value used in hydraulic design.
manning_n = 0.016

# Design rainfall intensity 'i' in mm/hr.
# Since IDF curves are not provided, a value typical for a 10-year, short-duration
# storm is assumed for designing drainage on an arterial road.
rainfall_intensity_mmhr = 150

# Step 2: Calculate derived parameters for the formula
# Calculate total flow path length (L) in meters
flow_path_length_m = num_lanes * lane_width_m

# Convert cross-slope (S) from percent to a decimal (m/m)
cross_slope_decimal = cross_slope_percent / 100

# Convert rainfall intensity (i) from mm/hr to m/s for use in the SI unit formula
# 1 mm/hr = 1 / (1000 * 3600) m/s
rainfall_intensity_ms = rainfall_intensity_mmhr / (1000 * 3600)

# Step 3: Apply the kinematic wave formula to calculate water film thickness
# The formula is: d = [ (n * i * L) / S^0.5 ]^0.6
# where d is depth in meters.

# Numerator of the inner term
numerator = manning_n * rainfall_intensity_ms * flow_path_length_m

# Denominator of the inner term
denominator = math.sqrt(cross_slope_decimal)

# Calculate water film thickness in meters
water_film_thickness_m = (numerator / denominator) ** 0.6

# Convert the final result to millimeters
water_film_thickness_mm = water_film_thickness_m * 1000

# Step 4: Print the results and the equation
print("--- Design Parameters ---")
print(f"Flow Path Length (L): {flow_path_length_m:.2f} m")
print(f"Cross-Slope (S): {cross_slope_decimal:.4f} m/m")
print(f"Manning's Roughness (n): {manning_n:.3f} (Assumed)")
print(f"Design Rainfall Intensity (i): {rainfall_intensity_mmhr} mm/hr (Assumed)")
print("\n--- Calculation ---")
print("The water film thickness (d) is calculated using the formula:")
print("d_mm = [ (n * i_ms * L) / (S^0.5) ]^0.6 * 1000\n")

print("Plugging in the numbers:")
# Display the equation with the values used
print(f"d_mm = [ ({manning_n} * {rainfall_intensity_ms:.6e} * {flow_path_length_m}) / ({cross_slope_decimal}**0.5) ]**0.6 * 1000")

print("\n--- Result ---")
print(f"The design water film thickness is: {water_film_thickness_mm:.2f} mm")

# Final answer in the specified format
final_answer = round(water_film_thickness_mm, 2)
# <<<4.14>>>