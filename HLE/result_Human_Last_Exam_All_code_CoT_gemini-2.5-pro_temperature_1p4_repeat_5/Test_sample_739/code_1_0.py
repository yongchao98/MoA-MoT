import math

# Step 1: Define road geometry parameters from the problem
num_lanes = 3
lane_width_m = 3.6
cross_slope_percent = 1.75

# Step 2: Calculate the total drainage path length (L)
drainage_path_length_m = num_lanes * lane_width_m

# Step 3: Assume a design rainfall intensity (i)
rainfall_intensity_mm_hr = 150.0

# Step 4 & 5: Apply the formula and calculate the result
# The formula is: WFT = 0.12 * L^0.36 * i^0.49 * S^-0.15
# where L is in meters, i is in mm/hr, and S is in percent.

# Calculate each term of the equation
term_L = math.pow(drainage_path_length_m, 0.36)
term_i = math.pow(rainfall_intensity_mm_hr, 0.49)
term_S = math.pow(cross_slope_percent, -0.15)

# Calculate the final water film thickness
water_film_thickness_mm = 0.12 * term_L * term_i * term_S

# Print the inputs and the final equation for clarity
print("--- Input Values ---")
print(f"Number of lanes: {num_lanes}")
print(f"Lane width: {lane_width_m} m")
print(f"Total drainage path length (L): {drainage_path_length_m:.2f} m")
print(f"Cross-slope (S): {cross_slope_percent}%")
print(f"Assumed rainfall intensity (i): {rainfall_intensity_mm_hr} mm/hr")
print("\n--- Calculation ---")
print("Using the formula: WFT = 0.12 * L^0.36 * i^0.49 * S^-0.15")
print(f"The equation with values is:")
print(f"WFT = 0.12 * ({drainage_path_length_m:.2f}^0.36) * ({rainfall_intensity_mm_hr}^0.49) * ({cross_slope_percent}^-0.15)")
print(f"WFT = 0.12 * {term_L:.4f} * {term_i:.4f} * {term_S:.4f}")
print("\n--- Result ---")
print(f"The design water film thickness is {water_film_thickness_mm:.2f} mm.")
<<<3.14>>>