import math

# Step 1: Define given parameters and constants
# Lane width in meters
lane_width = 3.6
# Number of lanes sloping to the edge
num_lanes = 3
# Cross-slope in m/m
cross_slope = 1.75 / 100
# Manning's roughness coefficient for rough-textured asphalt
n = 0.016

# Calculate the total width of the flow path (W_total)
W_total = num_lanes * lane_width

# Step 2: Assume a standard IDF curve for a 10-year storm
# The relationship is i = a / (T_c + b)^c, where T_c is in minutes and i is in mm/hr.
# These parameters are for a generic location and must be assumed.
idf_a = 998
idf_b = 6.05
idf_c = 0.814

# Constant for the Time of Concentration (Tc) formula in SI units
# Tc (min) = K * (n * W_total)^0.6 / (i^0.4 * S^0.3)
K = 0.938

print("This script calculates the design water film thickness based on an iterative solution.")
print("It assumes a 10-year storm IDF curve: i = 998 / (T_c + 6.05)^0.814.\n")

# Step 3: Iteratively solve for rainfall intensity (i)
# Start with an initial guess for Tc in minutes
T_c = 5.0

# Iterate 5 times to ensure convergence
for _ in range(5):
    # Calculate rainfall intensity 'i' (mm/hr) from the current Tc
    i = idf_a / ((T_c + idf_b) ** idf_c)
    
    # Calculate a new time of concentration 'T_c' (minutes) from the calculated 'i'
    numerator_tc = K * ((n * W_total) ** 0.6)
    denominator_tc = (i ** 0.4) * (cross_slope ** 0.3)
    T_c = numerator_tc / denominator_tc

# The final design rainfall intensity after convergence
design_i = idf_a / ((T_c + idf_b) ** idf_c)

# Step 4: Calculate the final water film thickness (d) using the converged intensity
# The formula is d(m) = [(n * i * W_total) / (3600000 * S^(1/2))]^(3/5)
numerator_d = n * design_i * W_total
denominator_d = 3600000 * (cross_slope ** 0.5)
d_meters = (numerator_d / denominator_d) ** (3.0/5.0)

# Convert the result to millimeters
d_mm = d_meters * 1000

# Step 5: Present the final result and the equation
print("--- Calculation Summary ---")
print(f"Converged Time of Concentration (Tc): {T_c:.2f} minutes")
print(f"Resulting Design Rainfall Intensity (i): {design_i:.2f} mm/hr\n")

print("--- Final Water Film Thickness Calculation ---")
print("The water film thickness (d) is calculated using the formula:")
print("d = [(n * i * W_total) / (3600000 * S^(1/2))]^(3/5)\n")
print("Plugging in the final values:")
print(f"d = [({n:.3f} * {design_i:.2f} mm/hr * {W_total:.1f} m) / (3600000 * {cross_slope:.4f}^(1/2))]^(3/5)")
print(f"d = [{numerator_d:.2f} / (3600000 * {math.sqrt(cross_slope):.5f})]^0.6")
print(f"d = [{numerator_d:.2f} / {denominator_d:.2f}]^0.6")
print(f"d = {d_meters:.6f} meters\n")

print("--- Final Answer ---")
print(f"The design water film thickness is {d_mm:.2f} mm.")

# The final answer in the requested format
final_answer = f"<<<{d_mm:.2f}>>>"