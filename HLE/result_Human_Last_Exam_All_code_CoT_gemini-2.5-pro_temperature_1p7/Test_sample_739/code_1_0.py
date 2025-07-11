import math

# Step 1: Define parameters from the problem statement and standard assumptions.

# Number of lanes in one direction
num_lanes = 3

# Width of each lane in meters
lane_width_m = 3.6

# Pavement cross-slope in percent
slope_percent = 1.75

# Manning's roughness coefficient for rough-textured asphalt (standard assumption)
n = 0.016

# Design rainfall intensity in mm/hr (standard assumption for pavement design)
I_mmhr = 150.0

# Step 2: Perform necessary calculations and unit conversions.

# Calculate the total drainage path length (L) in meters.
# Water flows across all three lanes from the median.
L_m = num_lanes * lane_width_m

# Convert the cross-slope from percent to decimal (m/m).
S_decimal = slope_percent / 100.0

# Convert rainfall intensity from mm/hr to m/s for the formula.
# 1 m = 1000 mm, 1 hr = 3600 s
i_ms = I_mmhr / (1000 * 3600)

# Step 3: Apply the water film thickness formula.
# d = [ (n * L * i) / (S^0.5) ]^0.6
# This calculates the depth in meters.

# Calculate the numerator of the fraction inside the brackets
numerator = n * L_m * i_ms

# Calculate the denominator
denominator = math.sqrt(S_decimal)

# Calculate the water film thickness in meters
d_m = (numerator / denominator) ** 0.6

# Convert the final result from meters to millimeters
d_mm = d_m * 1000

# Step 4: Print the results, including the equation with substituted values.
print("Calculation for Design Water Film Thickness")
print("-----------------------------------------")
print(f"Drainage Path Length (L): {L_m:.1f} m ({num_lanes} lanes @ {lane_width_m} m/lane)")
print(f"Pavement Cross-Slope (S): {S_decimal:.4f} m/m ({slope_percent}%)")
print(f"Manning's Roughness (n): {n} (for rough-textured asphalt)")
print(f"Design Rainfall Intensity (I): {I_mmhr} mm/hr -> {i_ms:.5e} m/s")
print("\nUsing the formula: d = [ (n * L * i) / (S^0.5) ]^0.6")
print("\nSubstituting the values into the formula:")
print(f"d_mm = [ ({n} * {L_m:.1f} * {i_ms:.5e}) / ({S_decimal:.4f}**0.5) ]**0.6 * 1000")

# Showing the intermediate calculation steps for clarity
print(f"d_mm = [ {numerator:.5e} / {denominator:.5f} ]**0.6 * 1000")
print(f"d_mm = [ {(numerator/denominator):.5e} ]**0.6 * 1000")
print(f"d_mm = {d_m:.5f} * 1000")

# Final Answer
print("\n-----------------------------------------")
print(f"The design water film thickness is: {d_mm:.2f} mm")
<<<3.88>>>