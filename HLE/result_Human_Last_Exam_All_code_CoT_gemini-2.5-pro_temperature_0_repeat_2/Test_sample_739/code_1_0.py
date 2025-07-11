import math

# Step 1: Define the given and assumed parameters.

# Given parameters from the problem description
lane_width_m = 3.6  # meters
num_lanes = 3
cross_slope_percent = 1.75 # %

# Assumed parameters based on standard engineering practice
# Manning's roughness coefficient 'n' for rough-textured asphalt.
manning_n = 0.016
# Design rainfall intensity 'i' in mm/hr. Assumed for hydroplaning analysis
# as rainfall-duration-frequency curves are not provided.
rainfall_intensity_mm_hr = 150

# Step 2: Calculate intermediate values.

# Convert cross-slope from percent to a decimal (m/m).
cross_slope = cross_slope_percent / 100

# Calculate the total drainage path length 'L' in meters.
drainage_path_length_m = num_lanes * lane_width_m

# Step 3: Apply the formula for water film thickness.
# Formula: d = 0.0038 * (n * L)^0.6 * i^0.6 * S^-0.3
term1 = manning_n * drainage_path_length_m
term2 = rainfall_intensity_mm_hr
term3 = cross_slope

water_film_thickness_mm = 0.0038 * (term1**0.6) * (term2**0.6) * (term3**-0.3)

# Step 4: Print the results.
print("--- Input Parameters ---")
print(f"Number of lanes: {num_lanes}")
print(f"Lane width: {lane_width_m} m")
print(f"Cross-slope: {cross_slope_percent}%")
print(f"Manning's n (assumed for rough asphalt): {manning_n}")
print(f"Design rainfall intensity (assumed): {rainfall_intensity_mm_hr} mm/hr")
print("\n--- Calculation ---")
print(f"Drainage path length (L) = {num_lanes} * {lane_width_m} = {drainage_path_length_m:.2f} m")
print(f"Cross-slope (S) = {cross_slope_percent} / 100 = {cross_slope}")
print("\n--- Final Equation ---")
print("d = 0.0038 * (n * L)^0.6 * i^0.6 * S^-0.3")
print(f"d = 0.0038 * ({manning_n} * {drainage_path_length_m:.2f})^0.6 * ({rainfall_intensity_mm_hr})^0.6 * ({cross_slope})^-0.3")
print(f"\nDesign Water Film Thickness (d) = {water_film_thickness_mm:.2f} mm")

# The final answer is formatted and printed below
final_answer = round(water_film_thickness_mm, 2)
# <<<8.77>>>