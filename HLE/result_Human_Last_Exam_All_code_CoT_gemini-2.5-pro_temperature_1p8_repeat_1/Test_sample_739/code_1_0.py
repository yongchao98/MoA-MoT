import math

# Step 1: Define input variables from the problem statement
lane_width_m = 3.6  # meters
num_lanes = 3
slope_percent = 1.75  # percent

# Step 2: State necessary engineering assumptions
# Manning's roughness coefficient for rough-textured asphalt
manning_n = 0.016
# Design rainfall intensity for hydroplaning analysis
rainfall_intensity_mm_hr = 100.0  # mm/hr

# Step 3: Convert inputs to SI base units for the formula
# Calculate the total drainage path length (L) in meters
drainage_length_L = num_lanes * lane_width_m
# Convert cross-slope (S) from percent to a decimal value (m/m)
slope_S = slope_percent / 100.0
# Convert rainfall intensity (I) from mm/hr to m/s
# 1 mm/hr = 1 / (1000 * 3600) m/s
rainfall_intensity_I_m_per_s = rainfall_intensity_mm_hr / (1000.0 * 3600.0)

# Step 4: Apply the Kinematic Wave equation
# T_d = ( (n * I * L) / (S^0.5) ) ^ 0.6
# This will calculate the thickness in meters
numerator = manning_n * rainfall_intensity_I_m_per_s * drainage_length_L
denominator = math.sqrt(slope_S)
water_film_thickness_m = (numerator / denominator)**0.6

# Step 5: Convert the final result to millimeters
water_film_thickness_mm = water_film_thickness_m * 1000.0

# Print the summary of parameters and the final calculation
print("--- Design Parameters ---")
print(f"Drainage Path Length (L): {drainage_length_L:.2f} m ({num_lanes} lanes @ {lane_width_m} m/lane)")
print(f"Pavement Cross-Slope (S): {slope_S:.4f} m/m ({slope_percent}%)")
print(f"Manning's Roughness (n): {manning_n:.3f} (Assumed for rough asphalt)")
print(f"Design Rainfall Intensity (I): {rainfall_intensity_mm_hr:.1f} mm/hr = {rainfall_intensity_I_m_per_s:.3e} m/s (Assumed for hydroplaning analysis)")
print("\n--- Calculation ---")
print("The water film thickness (T_d) is calculated using the Kinematic Wave equation:")
print("T_d (mm) = [ (n * I * L) / (S^0.5) ] ^ 0.6 * 1000")
print("\nPlugging in the values:")
print(f"T_d (mm) = [ ({manning_n:.3f} * {rainfall_intensity_I_m_per_s:.3e} * {drainage_length_L:.1f}) / ({slope_S:.4f}^0.5) ] ^ 0.6 * 1000")
print(f"T_d (mm) = [ {numerator:.3e} / {denominator:.4f} ] ^ 0.6 * 1000")
print(f"T_d (mm) = [ {(numerator / denominator):.3e} ] ^ 0.6 * 1000")
print(f"\nDesign Water Film Thickness = {water_film_thickness_mm:.2f} mm")
