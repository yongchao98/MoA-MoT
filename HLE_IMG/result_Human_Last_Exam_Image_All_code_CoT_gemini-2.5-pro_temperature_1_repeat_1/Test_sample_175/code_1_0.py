import math

# Step 1: Define the given parameters
B = 2.2  # Foundation width in m
L = 2.2  # Foundation length in m
D = 2.0  # Foundation depth in m
d_w = 0.6  # Groundwater level depth in m
gamma_dry = 18.0  # Unit weight of clay above groundwater in kN/m^3
gamma_sat = 20.0  # Unit weight of clay below groundwater in kN/m^3
c_uk = 57.0  # Characteristic undrained shear strength in kN/m^2 (kPa)
# For Design Approach 1, Combination 1 (DA1/1), the partial factor for undrained shear strength is 1.0
gamma_cu = 1.0

# Step 2: Calculate intermediate values
# Foundation Area (A')
A_prime = B * L

# Overburden pressure at foundation base (q)
# Pressure from soil above groundwater table + pressure from soil below groundwater table
q = (gamma_dry * d_w) + (gamma_sat * (D - d_w))

# Shape factor (s_c) for a square foundation in undrained analysis
s_c = 1.2

# Inclination factor (i_c) for vertical load
i_c = 1.0

# Design undrained shear strength (c_ud)
c_ud = c_uk / gamma_cu

# Step 3: Calculate the design bearing resistance (R_d)
# The formula for design bearing resistance in undrained conditions is:
# R_d = A' * [ (pi + 2) * c_ud * s_c * i_c + q ]
R_d = A_prime * ((math.pi + 2) * c_ud * s_c * i_c + q)

# Step 4: Print the results and the final equation
print("Calculation of Design Bearing Resistance (R_d) for Undrained Conditions")
print("-" * 70)
print(f"Given Parameters:")
print(f"  Foundation width B = {B} m")
print(f"  Foundation length L = {L} m")
print(f"  Foundation depth D = {D} m")
print(f"  Groundwater depth d_w = {d_w} m")
print(f"  Soil unit weight (above GWT) γ_dry = {gamma_dry} kN/m³")
print(f"  Soil unit weight (below GWT) γ_sat = {gamma_sat} kN/m³")
print(f"  Characteristic undrained shear strength c_uk = {c_uk} kN/m²")
print(f"  Partial factor for c_u (DA1/1) γ_cu = {gamma_cu}")
print("-" * 70)

print(f"Step A: Calculate Foundation Area (A')")
print(f"  A' = B * L = {B} * {L} = {A_prime:.2f} m²")
print("-" * 70)

print(f"Step B: Calculate Overburden Pressure (q)")
print(f"  q = (γ_dry * d_w) + (γ_sat * (D - d_w))")
print(f"  q = ({gamma_dry} * {d_w}) + ({gamma_sat} * ({D} - {d_w})) = {gamma_dry * d_w:.2f} + {gamma_sat * (D-d_w):.2f} = {q:.2f} kN/m²")
print("-" * 70)

print(f"Step C: Determine Factors and Design Parameters")
print(f"  Shape factor for square foundation s_c = {s_c}")
print(f"  Inclination factor for vertical load i_c = {i_c}")
print(f"  Design undrained shear strength c_ud = c_uk / γ_cu = {c_uk} / {gamma_cu} = {c_ud:.2f} kN/m²")
print("-" * 70)

print("Final Calculation of Design Resistance (R_d)")
print("R_d = A' * [ (π + 2) * c_ud * s_c * i_c + q ]")
# Printing the equation with all the values
print(f"R_d = {A_prime:.2f} * [ ({math.pi:.3f} + 2) * {c_ud:.1f} * {s_c:.1f} * {i_c:.1f} + {q:.1f} ]")
nc_term = (math.pi + 2) * c_ud * s_c * i_c
print(f"R_d = {A_prime:.2f} * [ {nc_term:.2f} + {q:.1f} ]")
print(f"R_d = {A_prime:.2f} * [ {nc_term + q:.2f} ]")
print(f"R_d = {R_d:.2f} kN")
print("\nThe design resistance under combination 1 for undrained condition is:")
print(f"{R_d:.1f} kN")
<<<1890.8>>>