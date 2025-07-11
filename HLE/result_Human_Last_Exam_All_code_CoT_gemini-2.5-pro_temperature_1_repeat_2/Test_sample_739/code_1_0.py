import math

# Step 1: Define given and assumed parameters
lanes = 3
lane_width_m = 3.6  # width of each lane in meters
slope_percent = 1.75  # cross-slope in percent
manning_n = 0.016  # Assumed Manning's roughness coefficient for rough asphalt
rainfall_intensity_mm_hr = 150  # Assumed design rainfall intensity in mm/hr

# Step 2: Calculate derived parameters and convert units
# Calculate the total flow path length (L) in meters
L_m = lanes * lane_width_m

# Convert the cross-slope (S) from percent to a decimal
S_decimal = slope_percent / 100

# Convert rainfall intensity (i) from mm/hr to m/s for the formula
# 1 mm/hr = 1 / (1000 * 3600) m/s
i_m_s = rainfall_intensity_mm_hr / (1000 * 3600)

# Step 3: Calculate the water film thickness using the formula
# The formula is d = [ (i * L * n) / (S^0.5) ]^0.6
# Calculate the numerator part: (i * L * n)
numerator = i_m_s * L_m * manning_n

# Calculate the denominator part: (S^0.5)
denominator = math.sqrt(S_decimal)

# Calculate the water film thickness in meters (d_m)
d_m = (numerator / denominator) ** 0.6

# Step 4: Convert the final result to millimeters (mm)
d_mm = d_m * 1000

# --- Output the results ---
print("--- Design Water Film Thickness Calculation ---")
print("\nParameters Used:")
print(f"  - Number of lanes: {lanes}")
print(f"  - Lane width: {lane_width_m} m")
print(f"  - Total flow path length (L): {L_m:.1f} m")
print(f"  - Cross-slope (S): {slope_percent}% or {S_decimal:.4f}")
print(f"  - Manning's n (assumed): {manning_n}")
print(f"  - Design rainfall intensity (i, assumed): {rainfall_intensity_mm_hr} mm/hr or {i_m_s:.3e} m/s")

print("\nFormula:")
print("  d = [ (i * L * n) / (S^0.5) ]^0.6 * 1000")

print("\nCalculation with values:")
# Using f-string formatting to display the numbers in the equation
print(f"  d = [ ({i_m_s:.6f} * {L_m:.1f} * {manning_n:.3f}) / ({S_decimal:.4f}^0.5) ]^0.6 * 1000")
print(f"  d = [ ({numerator:.7f}) / ({denominator:.5f}) ]^0.6 * 1000")
print(f"  d = [ ({(numerator/denominator):.7f}) ]^0.6 * 1000")
print(f"  d = {d_m:.6f} * 1000")

print("\n--- Result ---")
print(f"The calculated design water film thickness is: {d_mm:.2f} mm")

# Final answer in the specified format
print(f"\n<<<{d_mm:.2f}>>>")