import math

# Step 1: Define the given and assumed parameters
# Given parameters from the problem description
num_lanes = 3
lane_width_m = 3.6  # Lane width in meters
cross_slope_percent = 1.75  # Cross-slope in percent

# Assumed parameters for hydroplaning analysis
# Manning's roughness coefficient 'n' for rough-textured asphalt
n_roughness = 0.016
# Design rainfall intensity 'i' in mm/hr. A high intensity is chosen for critical hydroplaning checks.
i_mm_per_hr = 150.0

# Step 2: Calculate variables for the formula
# The total flow path length 'L' is the width of the three lanes
L_m = num_lanes * lane_width_m

# The cross-slope 'S' as a decimal
S_slope = cross_slope_percent / 100.0

# Convert rainfall intensity 'i' from mm/hr to m/s for use in the SI unit-based formula
i_m_per_s = i_mm_per_hr / (1000 * 3600)

# Step 3: Calculate the water film thickness using the formula
# The formula is: d = [(i * L * n) / S^(1/2)]^(3/5)
numerator = i_m_per_s * L_m * n_roughness
denominator = math.sqrt(S_slope)
term_inside_brackets = numerator / denominator

# Water film thickness in meters
d_m = term_inside_brackets**(3.0/5.0)

# Convert the final result to millimeters
d_mm = d_m * 1000

# Step 4: Print the process and the final result
print("Calculation of Design Water Film Thickness\n")
print("The formula used is: d = [(i * L * n) / S^(1/2)]^(3/5)\n")

print("Parameters used in the calculation:")
print(f"  - Flow Path Length (L): {L_m:.2f} m ({num_lanes} lanes * {lane_width_m} m/lane)")
print(f"  - Cross-Slope (S): {cross_slope_percent}% = {S_slope}")
print(f"  - Manning's Roughness (n): {n_roughness} (Assumed for rough-textured asphalt)")
print(f"  - Design Rainfall Intensity (i): {i_mm_per_hr} mm/hr = {i_m_per_s:.8f} m/s (Assumed for hydroplaning analysis)\n")

print("Substituting the values into the formula:")
# The user request: "output each number in the final equation"
print(f"d (m) = [({i_m_per_s:.8f} * {L_m:.2f} * {n_roughness}) / ({S_slope})^(1/2)]^(3/5)")
print(f"d (m) = [({numerator:.8f}) / ({denominator:.5f})]^(3/5)")
print(f"d (m) = ({term_inside_brackets:.8f})^(3/5)")
print(f"d (m) = {d_m:.5f} m\n")

print("Converting the result to millimeters:")
print(f"d (mm) = {d_m:.5f} m * 1000 mm/m\n")

print(f"The final design water film thickness is {d_mm:.2f} mm.")

<<<5.57>>>