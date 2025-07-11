import math

# Step 1: Define given parameters and assumptions
num_lanes = 3
lane_width_m = 3.6
cross_slope_percent = 1.75

# Assumptions for parameters not explicitly provided
# Manning's roughness coefficient (n) for rough-textured asphalt
manning_n = 0.016
# Design rainfall intensity (i) in mm/hr (a common value for hydroplaning analysis)
rainfall_intensity_mm_hr = 150

print("--- Parameters and Assumptions ---")
print(f"Number of lanes per direction: {num_lanes}")
print(f"Lane width: {lane_width_m} m")
print(f"Pavement cross-slope: {cross_slope_percent}%")
print(f"Assumed Manning's roughness coefficient (n) for rough asphalt: {manning_n}")
print(f"Assumed design rainfall intensity (i): {rainfall_intensity_mm_hr} mm/hr")
print("-" * 34 + "\n")


# Step 2: Prepare variables for the formula
# Calculate the total drainage path length (L) across the three lanes
L_m = num_lanes * lane_width_m

# Convert cross-slope from percent to decimal
S = cross_slope_percent / 100.0

# Convert rainfall intensity from mm/hr to m/s for SI unit consistency
# 1 mm/hr = 1/1000 m/hr = 1/(1000 * 3600) m/s
i_m_s = rainfall_intensity_mm_hr / (1000 * 3600)


# Step 3: Apply the Kinematic Wave Equation
# d = (n * L)^(3/5) * i^(3/5) * S^(-3/10)
# Calculate each component of the formula
term1 = (manning_n * L_m)**(3/5)
term2 = i_m_s**(3/5)
term3 = S**(-3/10)

# Calculate the water film thickness in meters
d_m = term1 * term2 * term3

# Convert the final result to millimeters
d_mm = d_m * 1000


# Step 4: Display the results
print("--- Calculation ---")
print("The water film thickness (d) is calculated using the Kinematic Wave Equation:")
print("d = (n * L)^(3/5) * i^(3/5) * S^(-3/10)\n")
print("Substituting the values into the equation:")
# We display the numerical equation as requested
print(f"d = ({manning_n} * {L_m:.1f})^(3/5) * ({i_m_s:.2e})^(3/5) * ({S:.4f})^(-3/10)")
print(f"d = ({term1:.5f}) * ({term2:.5f}) * ({term3:.5f})")
print(f"d = {d_m:.5f} m\n")

print("--- Final Answer ---")
print(f"The calculated design water film thickness is: {d_mm:.3f} mm")

# <<<1.232>>>