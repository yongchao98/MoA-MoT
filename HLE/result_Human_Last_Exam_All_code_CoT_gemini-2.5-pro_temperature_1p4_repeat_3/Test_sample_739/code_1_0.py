import math

# --- Given Parameters ---
num_lanes = 3
lane_width_m = 3.6  # meters
cross_slope_percent = 1.75

# --- Assumed Parameters (based on standard engineering practice) ---
# Manning's roughness coefficient 'n' for rough-textured asphalt
manning_n = 0.016
# Design rainfall intensity 'i' in mm/hr.
# This is a critical assumption as IDF curves were not provided.
# A value of 150 mm/hr represents a heavy design storm suitable for hydroplaning analysis.
rainfall_intensity_mm_hr = 150.0

# --- Step-by-step Calculation ---

# 1. Calculate the total flow path length (L) in meters
# This is the total width of the travelled way in one direction.
L_m = num_lanes * lane_width_m

# 2. Convert parameters to consistent SI units for the formula
# Convert cross-slope from percent to a decimal (m/m)
S = cross_slope_percent / 100.0

# Convert rainfall intensity from mm/hr to meters/second
i_m_per_s = rainfall_intensity_mm_hr / (1000 * 3600)

# 3. Apply the kinematic wave equation to find water film thickness in meters
# d = (i * L * n)^0.6 * S^(-0.3)
# Note: This is equivalent to d = ((i * L * n) / (S^0.5))^0.6
term1 = i_m_per_s * L_m * manning_n
term2 = S**(-0.3)
d_m = (term1**0.6) * term2

# 4. Convert the final result to millimeters
d_mm = d_m * 1000

# --- Output the results ---
print("--- Input and Assumed Parameters ---")
print(f"Number of lanes: {num_lanes}")
print(f"Lane width: {lane_width_m} m")
print(f"Total flow path length (L): {L_m:.2f} m")
print(f"Cross-slope (S): {cross_slope_percent}% or {S:.4f} m/m")
print(f"Assumed Manning's n: {manning_n}")
print(f"Assumed Rainfall Intensity (i): {rainfall_intensity_mm_hr} mm/hr or {i_m_per_s:.8f} m/s")
print("\n--- Calculation ---")
print("The formula for water film thickness (d) in meters is: d = (i * L * n)^0.6 * S^(-0.3)")
print("Substituting the values in SI units (meters and seconds):")
equation_str = (
    f"d_mm = (({i_m_per_s:.8f} * {L_m:.2f} * {manning_n})**0.6 * "
    f"{S:.4f}**(-0.3)) * 1000"
)
print(equation_str)
print(f"\nCalculated water film thickness: {d_mm:.2f} mm")
