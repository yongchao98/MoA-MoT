import math

# Step 1: Define the input parameters from the problem description.
num_lanes = 3
lane_width_m = 3.6
cross_slope_percent = 1.75
# For rough-textured asphalt, a standard Manning's n value is 0.016.
manning_n = 0.016
# A standard design rainfall intensity for hydroplaning analysis is 150 mm/hr.
rainfall_intensity_mm_hr = 150

# Step 2: Calculate derived parameters and convert units for the formula.
# Calculate the total flow path length (L) in meters.
L_m = num_lanes * lane_width_m
# Convert the cross-slope (S) from percent to a decimal.
S_decimal = cross_slope_percent / 100
# Convert rainfall intensity (i) from mm/hr to m/s for use in the metric formula.
i_m_s = rainfall_intensity_mm_hr / (1000 * 3600)

# Step 3: Apply the kinematic wave formula to calculate water film thickness.
# The formula is: d = [(n * i * L) / S^0.5]^(3/5) or d = [(n * i * L) / S^0.5]^0.6
# All variables are in base SI units (meters, seconds).
try:
    # Numerator of the term inside the brackets
    numerator = manning_n * i_m_s * L_m
    # Denominator of the term inside the brackets
    denominator = math.sqrt(S_decimal)
    
    # Calculate water film thickness in meters
    d_m = (numerator / denominator) ** 0.6
    
    # Convert the final result to millimeters
    d_mm = d_m * 1000
    
    # Step 4: Print the equation with values and the final result.
    print("The design water film thickness (d) is calculated using the kinematic wave equation:")
    print("d (m) = [(n * i * L) / (S^0.5)]^0.6\n")
    print("Where:")
    print(f"  n (Manning's Roughness) = {manning_n}")
    print(f"  i (Rainfall Intensity) = {rainfall_intensity_mm_hr} mm/hr = {i_m_s:.6f} m/s")
    print(f"  L (Flow Path Length) = {L_m} m")
    print(f"  S (Cross-Slope) = {cross_slope_percent}% = {S_decimal}\n")
    
    print("Substituting the values into the equation:")
    print(f"d (m) = [({manning_n} * {i_m_s:.6f} * {L_m}) / ({S_decimal}^0.5)]^0.6")
    print(f"d (m) = [({numerator:.6e}) / ({denominator:.4f})]^0.6")
    print(f"d (m) = [{(numerator / denominator):.6e}]^0.6")
    print(f"d (m) = {d_m:.6f}\n")
    
    print("Final Answer:")
    print(f"The design water film thickness is {d_mm:.2f} mm.")

except (ValueError, ZeroDivisionError) as e:
    print(f"An error occurred during calculation: {e}")
