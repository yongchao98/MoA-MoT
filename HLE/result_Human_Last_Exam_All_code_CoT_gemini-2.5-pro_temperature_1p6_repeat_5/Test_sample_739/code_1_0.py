import math

# --- 1. Define Parameters ---
# Given parameters
lane_width_m = 3.6  # m
num_lanes = 3
slope_percent = 1.75  # %

# Assumed parameters based on engineering standards
# Manning's roughness coefficient (n) for "rough-textured asphalt pavement".
# A standard value for dense-graded asphalt pavement is used.
n = 0.016

# Design rainfall intensity (i). Since IDF curves are not provided, a typical
# high-intensity value for arterial road drainage design is assumed.
i_mm_hr = 150  # mm/hr

# --- 2. Preparatory Calculations ---
# Calculate total drainage path length (L) from the center to the outer edge
L = num_lanes * lane_width_m

# Convert slope (S) from percent to decimal
S = slope_percent / 100

# Convert rainfall intensity (i) from mm/hr to m/s for use in the formula
i_m_s = i_mm_hr / (1000 * 3600)  # (mm/hr) * (1m/1000mm) * (1hr/3600s)

# --- 3. Main Calculation using the Kinematic Wave Formula ---
# The formula for water film thickness (dw) is: dw = [ (n * L * i) / S^0.5 ] ^ 0.6

# Calculate the numerator and denominator of the base term
numerator = n * L * i_m_s
denominator = math.sqrt(S)
base = numerator / denominator

# Calculate the water film thickness in meters
dw_m = base ** 0.6

# Convert the final result to millimeters
dw_mm = dw_m * 1000

# --- 4. Output the Results ---
print("--- Design Parameters ---")
print(f"Total drainage path length (L): {L:.1f} m")
print(f"Pavement cross-slope (S): {S:.4f} m/m")
print(f"Assumed Manning's roughness coefficient (n): {n}")
print(f"Assumed design rainfall intensity (i): {i_mm_hr} mm/hr")
print("")

print("--- Calculation Steps ---")
print("The water film thickness (dw) is calculated using the formula: dw = [ (n * L * i) / (S^0.5) ] ^ 0.6")
print(f"\n1. Substitute values into the equation (with i in m/s = {i_m_s:.8f}):")
print(f"dw_meters = [ ({n} * {L:.1f} * {i_m_s:.8f}) / ({S:.4f}^0.5) ] ^ 0.6")
print("")

print(f"2. Calculate the terms inside the brackets:")
print(f"dw_meters = [ {numerator:.8f} / {denominator:.6f} ] ^ 0.6")
print(f"dw_meters = [ {base:.8f} ] ^ 0.6")
print("")

print(f"3. Perform the exponentiation to find the depth in meters:")
print(f"dw_meters = {dw_m:.6f} m")
print("")

print("--- Final Answer ---")
print(f"4. Convert the water film thickness to millimeters:")
print(f"dw_mm = {dw_m:.6f} m * 1000 mm/m")
print(f"The design water film thickness is {dw_mm:.2f} mm.")

<<<4.16>>>