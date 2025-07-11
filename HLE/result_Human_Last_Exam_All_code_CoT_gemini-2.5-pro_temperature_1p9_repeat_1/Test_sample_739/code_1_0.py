import math

# --- Step 1: Define parameters from the problem and state assumptions ---

# Pavement and geometric properties
lane_width_m = 3.6
num_lanes = 3
cross_slope_percent = 1.75

# Assumed parameters based on standard engineering practice
# Manning's roughness coefficient 'n' for rough-textured asphalt
n_manning = 0.015
# Design rainfall intensity 'i' in mm/hr. Assumed because IDF curves were not provided.
i_mm_hr = 125

# --- Step 2: Calculate intermediate values and convert units ---

# L: Total length of the drainage path across the traveled way in meters
L_path_m = num_lanes * lane_width_m

# S: Cross-slope as a decimal (m/m)
S_slope = cross_slope_percent / 100

# i: Convert rainfall intensity from mm/hr to m/s for the formula
# 1 mm/hr = 1 / (1000 * 3600) m/s
i_m_s = i_mm_hr / (1000 * 3600)

# --- Step 3: Apply the Kinematic Wave Equation ---

# The formula for water film thickness (Td) in meters is:
# Td = [ (n * L * i) / (S^0.5) ]^(3/5)

# Calculate the terms for the equation
numerator = n_manning * L_path_m * i_m_s
sqrt_S = math.sqrt(S_slope)
base_of_exponent = numerator / sqrt_S
exponent = 3.0 / 5.0
Td_meters = math.pow(base_of_exponent, exponent)

# Convert the final result from meters to millimeters
Td_mm = Td_meters * 1000

# --- Step 4: Output the calculation steps and final answer ---

print("--- Water Film Thickness Calculation ---")
print(f"\nParameters:")
print(f"Flow Path Length (L) = {num_lanes} lanes * {lane_width_m} m/lane = {L_path_m:.1f} m")
print(f"Cross-Slope (S) = {cross_slope_percent}% = {S_slope}")
print(f"Assumed Manning's n = {n_manning}")
print(f"Assumed Rainfall Intensity (i) = {i_mm_hr} mm/hr = {i_m_s:.5e} m/s")

print("\nCalculation using Td = [ (n * L * i) / (S^0.5) ]^0.6:")
print(f"Td (m) = [ ({n_manning} * {L_path_m:.1f} * {i_m_s:.5e}) / ({S_slope}^0.5) ]^0.6")
print(f"Td (m) = [ ({numerator:.5e}) / ({sqrt_S:.5f}) ]^0.6")
print(f"Td (m) = [ {base_of_exponent:.5e} ]^0.6")
print(f"Td (m) = {Td_meters:.5f} m")

print("\n--- Final Answer ---")
print(f"The design water film thickness is {Td_mm:.2f} mm.")
<<<3.5>>>