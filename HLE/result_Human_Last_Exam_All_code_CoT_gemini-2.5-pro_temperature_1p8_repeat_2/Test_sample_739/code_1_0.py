import math

# Step 1: Define the parameters based on the problem statement and engineering standards.

# Number of lanes in one direction
num_lanes = 3
# Width of a single lane in meters
lane_width_m = 3.6
# Cross-slope of the pavement in percent
cross_slope_percent = 1.75
# Manning's roughness coefficient for rough-textured asphalt (standard value)
manning_n = 0.016
# Design rainfall intensity in mm/hr (standard assumption for hydroplaning analysis)
rainfall_intensity_mmhr = 100

# Step 2: Calculate the parameters for the formula in consistent SI units.

# L: Total drainage path length in meters
L = num_lanes * lane_width_m

# S: Cross-slope as a decimal
S = cross_slope_percent / 100

# i: Rainfall intensity in meters per second
# Conversion: 1 mm/hr = 1/(1000 * 3600) m/s
i_ms = rainfall_intensity_mmhr / (1000 * 3600)

# Step 3: Apply the kinematic wave equation to find the water film thickness.
# Formula: d_w = [ (n * L * i) / S^(1/2) ]^(3/5)

# Calculate the numerator and denominator separately for clarity
numerator = manning_n * L * i_ms
denominator = math.sqrt(S)
base = numerator / denominator

# Calculate the water film thickness in meters
d_w_meters = base ** (3/5)

# Convert the result to millimeters
d_w_mm = d_w_meters * 1000

# Step 4: Print the final equation and the result.
print("Calculation of Water Film Thickness (d_w)\n")
print("Formula: d_w (mm) = [ (n * L * i) / S^(1/2) ]^(3/5) * 1000\n")
print("Substituting the values:")
# We use the full precision values for calculation but show rounded for display
# Print the equation with values plugged in
equation_str = (
    f"d_w (mm) = [ ({manning_n} * {L:.1f} m * {i_ms:.8f} m/s) / "
    f"({S:.4f}^(1/2)) ]^(3/5) * 1000"
)
print(equation_str)
print("\nCalculating the final result:")
print(f"The design water film thickness is: {d_w_mm:.2f} mm")
<<<3.98>>>