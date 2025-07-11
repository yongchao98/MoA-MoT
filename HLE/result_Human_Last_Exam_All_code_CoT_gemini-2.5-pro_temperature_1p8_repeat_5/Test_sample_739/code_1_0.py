import math

# --- Step 1: Define Parameters ---
# The problem provides the road geometry and surface type. We must define these
# as variables and select standard values for unspecified engineering parameters.

# The road has three lanes in the direction of travel.
num_lanes = 3
# Each lane is 3.6 meters wide.
lane_width_m = 3.6
# The cross-slope is 1.75%, sloping outward.
cross_slope_percent = 1.75
# The pavement is rough-textured asphalt. For this surface, a standard
# Manning's roughness coefficient 'n' for sheet flow analysis is 0.016.
manning_n = 0.016
# Rainfall intensity 'I' must be assumed. For hydroplaning analysis, a high-intensity,
# short-duration storm (e.g., 10-year, 5-10 minute duration) is appropriate.
# A value of 150 mm/hr is a common and conservative assumption for this type of design.
rainfall_intensity_mmhr = 150.0

# --- Step 2: Calculate Derived Values and State the Formula ---
# Calculate the total flow path length (L) across the lanes.
L_m = num_lanes * lane_width_m
# Convert the cross-slope from a percentage to a decimal value (m/m).
S = cross_slope_percent / 100.0

# We will use the Papadakis-Wylie version of the kinematic wave equation, converted to SI units.
# The formula is: d_w = C * (n * L)^0.6 * I^0.4 * S^-0.2
# where d_w is the water film thickness in mm, and C is a unit conversion constant of ~0.053.
C = 0.05297

# --- Step 3: Perform the Calculation ---
# Calculate each term of the equation.
term1_base = manning_n * L_m
term2_base = rainfall_intensity_mmhr
term3_base = S

power_term1 = math.pow(term1_base, 0.6)
power_term2 = math.pow(term2_base, 0.4)
power_term3 = math.pow(term3_base, -0.2)

# Calculate the final water film thickness.
water_film_thickness_mm = C * power_term1 * power_term2 * power_term3

# --- Step 4: Output the Result ---
# The prompt requires printing the equation with the numbers substituted in.

print("--- Design Water Film Thickness Calculation ---")
print("\n[Input Values]")
print(f"Number of lanes: {num_lanes}")
print(f"Flow Path Length (L): {num_lanes} lanes * {lane_width_m} m/lane = {L_m:.2f} m")
print(f"Cross-Slope (S): {cross_slope_percent}% = {S:.4f} m/m")
print(f"Manning's Roughness (n): {manning_n} (for rough asphalt)")
print(f"Assumed Rainfall Intensity (I): {rainfall_intensity_mmhr} mm/hr")

print("\n[Formula]")
print("d_w = C * (n * L)^0.6 * I^0.4 * S^-0.2")

print("\n[Calculation Steps]")
print(f"d_w = 0.053 * ({manning_n} * {L_m:.2f})^0.6 * ({rainfall_intensity_mmhr})^0.4 * ({S})^-0.2")
print(f"d_w = 0.053 * ({term1_base:.4f})^0.6 * ({term2_base})^0.4 * ({term3_base})^-0.2")
print(f"d_w = 0.053 * {power_term1:.4f} * {power_term2:.4f} * {power_term3:.4f}")

print("\n[Final Result]")
print(f"The design water film thickness is {water_film_thickness_mm:.2f} mm.")
