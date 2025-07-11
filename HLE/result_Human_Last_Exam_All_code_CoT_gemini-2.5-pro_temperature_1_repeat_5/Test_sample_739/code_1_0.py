import math

# Step 1: Define the parameters based on the problem description and standard values.

# Manning's roughness coefficient for rough-textured asphalt
n = 0.016

# Flow path length (3 lanes * 3.6 m/lane)
num_lanes = 3
lane_width_m = 3.6
L = num_lanes * lane_width_m

# Pavement cross-slope (1.75%)
cross_slope_percent = 1.75
S = cross_slope_percent / 100

# Assumed design rainfall intensity in mm/hr (as no IDF curve was provided)
# A value of 125 mm/hr is a reasonable assumption for hydroplaning analysis.
i_mm_per_hr = 125

# Convert rainfall intensity to SI units (m/s)
# 1 hr = 3600 s, 1 m = 1000 mm
i_m_per_s = i_mm_per_hr / (1000 * 3600)

# Step 2: Calculate the water film thickness using the Kinematic Wave equation.
# d = [ (n * L * i) / (S^0.5) ]^(3/5)

# Calculate the term inside the brackets
val_inside_bracket = (n * L * i_m_per_s) / (S**0.5)

# Calculate the final water film thickness in meters
d_m = val_inside_bracket**0.6

# Convert the result to millimeters
d_mm = d_m * 1000

# Step 3: Print the final equation with all the numbers and the result.
print("Calculation of Water Film Thickness (d)")
print("=========================================")
print(f"Formula: d = [ (n * L * i) / (S^0.5) ]^0.6")
print("\nSubstituting the values:")
print(f"n (Manning's n) = {n}")
print(f"L (Flow Path Length) = {L} m")
print(f"i (Rainfall Intensity) = {i_mm_per_hr} mm/hr = {i_m_per_s:.8f} m/s")
print(f"S (Cross-slope) = {S}")
print("\nFinal Equation:")
# The prompt requires printing each number in the final equation.
final_equation = f"d = [ ({n} * {L} * {i_m_per_s:.8f}) / ({S}**0.5) ]^0.6"
print(final_equation)
print(f"d = [ {n * L * i_m_per_s:.8f} / {S**0.5:.6f} ]^0.6")
print(f"d = [ {val_inside_bracket:.8f} ]^0.6")
print(f"d = {d_m:.6f} m")
print("\n---")
print(f"Design Water Film Thickness = {d_mm:.2f} mm")
print("---")
<<<4.65>>>