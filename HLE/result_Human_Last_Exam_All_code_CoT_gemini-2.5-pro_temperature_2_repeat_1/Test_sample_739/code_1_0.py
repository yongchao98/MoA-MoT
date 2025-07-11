import math

# Step 1: Define the input parameters
lane_width = 3.6  # meters
num_lanes = 3
cross_slope_percent = 1.75  # %
# Manning's roughness coefficient for rough-textured asphalt
n = 0.016
# Assumed minimum time of concentration for design (duration of rainfall)
Tc_min = 5  # minutes

# Step 2: Calculate derived geometric parameters
# Flow path length (L) is the total width of the traveled way
L = num_lanes * lane_width
# Cross-slope (S) as a decimal
S = cross_slope_percent / 100.0

# Step 3: Determine the design rainfall intensity (I)
# Since IDF curves are not provided, we use a representative formula for a 10-year storm.
# Formula: I(mm/hr) = a / (Tc + b)^c. We use a=2500, b=15, c=0.85 as generic values.
I_mm_per_hr = 2500 / ((Tc_min + 15) ** 0.85)
# Convert rainfall intensity to m/s for use in the kinematic wave equation
# 1 mm/hr = 1 / (1000 * 3600) m/s
I_m_per_s = I_mm_per_hr / (1000 * 3600)

# Step 4: Calculate the water film thickness (d_w) using the Kinematic Wave equation
# d_w = [ (n * I * L) / S^0.5 ]^(3/5)
numerator = n * I_m_per_s * L
denominator = math.sqrt(S)
base = numerator / denominator
exponent = 3.0 / 5.0
d_w_m = base ** exponent

# Convert final answer to millimeters
d_w_mm = d_w_m * 1000

# Step 5: Print the results and the final equation
print("--- Design Parameters ---")
print(f"Flow Path Length (L): {L:.2f} m")
print(f"Cross-Slope (S): {S:.4f} m/m")
print(f"Manning's n: {n:.3f}")
print(f"Design Rainfall Duration (Tc): {Tc_min} min")
print(f"Design Rainfall Intensity (I): {I_mm_per_hr:.2f} mm/hr ({I_m_per_s:.6f} m/s)")
print("\n--- Water Film Thickness Calculation ---")
print("Equation: d_w = [ (n * I * L) / S^0.5 ]^(3/5)")
print(f"Values: d_w = [ ({n:.3f} * {I_m_per_s:.6f} * {L:.2f}) / {S:.4f}^0.5 ]^(3/5)")
print(f"Result: d_w = [ {numerator:.8f} / {denominator:.5f} ]^({exponent:.1f})")
print(f"Result: d_w = [ {base:.8f} ]^({exponent:.1f})")
print(f"\nThe design water film thickness is {d_w_m:.5f} m.")
print(f"Final Answer: The design water film thickness is {d_w_mm:.2f} mm.")

<<<2.89>>>