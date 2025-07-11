import math

# Step 1: Define constants and parameters from the problem statement
num_lanes = 3
lane_width_m = 3.6
cross_slope_percent = 1.75

# Step 2: Select a standard Manning's n for the pavement type
manning_n = 0.016  # Standard value for rough-textured asphalt

# Derived parameters
drainage_path_length_L = num_lanes * lane_width_m
slope_S = cross_slope_percent / 100

# Step 3: Define the assumed rainfall Intensity-Duration-Frequency (IDF) relationship.
# Since no IDF curve was provided, a typical relationship for a 10-year storm is assumed:
# i(mm/hr) = 3000 / (Tc(min) + 20)
def get_rainfall_intensity(tc_minutes):
    """Calculates rainfall intensity based on time of concentration."""
    return 3000.0 / (tc_minutes + 20.0)

# Step 4: Perform the iterative calculation
# Initialize Time of Concentration (Tc) with a reasonable guess
tc_minutes = 5.0
num_iterations = 10  # Number of iterations is sufficient for convergence

# Final values will be stored in these variables
final_rainfall_intensity_mmhr = 0
water_depth_m = 0

for i in range(num_iterations):
    # a) Calculate rainfall intensity 'i' (mm/hr) from Tc
    final_rainfall_intensity_mmhr = get_rainfall_intensity(tc_minutes)

    # b) Calculate water film depth 'dw' (m) using the Kinematic Wave equation
    # dw = ( (n * (i / 3600000) * L) / S^0.5 )^0.6
    i_m_per_sec = final_rainfall_intensity_mmhr / 3600000.0
    slope_sqrt = math.sqrt(slope_S)
    
    numerator = manning_n * i_m_per_sec * drainage_path_length_L
    water_depth_m = (numerator / slope_sqrt) ** 0.6

    # c) Calculate flow velocity 'v' (m/s) using Manning's equation for shallow flow
    # v = (1/n) * R^(2/3) * S^(1/2), where hydraulic radius R is approx. water_depth_m
    velocity_m_per_s = (1 / manning_n) * (water_depth_m ** (2.0 / 3.0)) * slope_sqrt
    
    # d) Update Tc based on the new flow velocity
    # If velocity is zero, avoid division by zero
    if velocity_m_per_s > 0:
      tc_seconds = drainage_path_length_L / velocity_m_per_s
      tc_minutes = tc_seconds / 60.0

# Convert final water depth to mm for the answer
water_depth_mm = water_depth_m * 1000

# Step 5: Print the detailed results
print("--- Design Water Film Thickness Calculation ---")

print("\n1. Input Parameters:")
print(f"  - Number of lanes per direction: {num_lanes}")
print(f"  - Lane width: {lane_width_m:.1f} m")
print(f"  - Total drainage path length (L): {drainage_path_length_L:.2f} m")
print(f"  - Cross-slope (S): {cross_slope_percent:.2f}% or {slope_S:.4f} m/m")
print(f"  - Manning's roughness 'n' (rough asphalt): {manning_n:.3f}")

print("\n2. Assumed Rainfall Data (as curves were not provided):")
print(f"  - IDF Relation: i(mm/hr) = 3000 / (Tc(min) + 20)")

print("\n3. Converged Intermediate Results:")
print(f"  - Time of Concentration (Tc): {tc_minutes:.2f} minutes")
print(f"  - Design Rainfall Intensity (i): {final_rainfall_intensity_mmhr:.2f} mm/hr")

print("\n4. Final Water Film Thickness (dw) Calculation:")
print(f"  - Using the Kinematic Wave Equation: dw(m) = ((n * i_mps * L) / S^0.5)^0.6")
print(f"  - With final values: dw(m) = (({manning_n:.3f} * ({final_rainfall_intensity_mmhr:.2f} / 3600000) * {drainage_path_length_L:.2f}) / {slope_S:.4f}^0.5)^0.6")

print("\n--- FINAL ANSWER ---")
print(f"The calculated design water film thickness is {water_depth_mm:.2f} mm.")
