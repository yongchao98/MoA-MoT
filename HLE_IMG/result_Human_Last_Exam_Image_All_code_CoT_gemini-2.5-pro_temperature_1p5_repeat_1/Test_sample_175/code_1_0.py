import math

# --- 1. Define given parameters ---
# Foundation dimensions
B = 2.2  # m
L = 2.2  # m, it's a square foundation

# Foundation and soil profile depths
D_f = 2.0  # m, depth to foundation base
d_w = 0.6  # m, depth to groundwater level

# Soil properties
c_uk = 57.0  # kN/m^2, characteristic undrained shear strength
gamma_soil_above = 18.0  # kN/m^3, unit weight of clay above GWL
gamma_soil_below = 20.0  # kN/m^3, unit weight of clay below GWL

# Partial factor of safety from Eurocode 7 (DA1-C1) for undrained bearing resistance
gamma_R = 1.4

# --- 2. Calculate intermediate values ---
# Foundation Area (A)
A = B * L

# Total overburden pressure at foundation base (q)
# Pressure from soil above groundwater level
q_part1 = d_w * gamma_soil_above
# Pressure from soil below groundwater level
q_part2 = (D_f - d_w) * gamma_soil_below
q = q_part1 + q_part2

# Factors for bearing capacity equation (undrained conditions)
# Bearing capacity factor Nc is (pi + 2)
Nc = math.pi + 2
# Shape factor (s_c) for a square foundation
s_c = 1.2
# Inclination factor (i_c) for vertical loading
i_c = 1.0

# --- 3. Calculate Characteristic Gross Bearing Resistance (R_k) ---
# Gross bearing pressure
q_gross_k = Nc * c_uk * s_c * i_c + q
# Total resistance
R_k = A * q_gross_k

# --- 4. Calculate Design Bearing Resistance (R_d) ---
R_d = R_k / gamma_R

# --- 5. Print the results and the equation ---
print("Calculation of the Design Bearing Resistance (R_d)")
print("--------------------------------------------------")
print("The formula for design resistance R_d is:")
print("R_d = R_k / gamma_R")
print("Where R_k = A * ((pi + 2) * c_uk * s_c * i_c + q)\n")

print("Values used in the calculation:")
print(f"Foundation Area, A = {B:.2f} m * {L:.2f} m = {A:.2f} m^2")
print(f"Overburden pressure, q = ({d_w:.2f} m * {gamma_soil_above:.2f} kN/m^3) + (({D_f:.2f} - {d_w:.2f}) m * {gamma_soil_below:.2f} kN/m^3) = {q:.2f} kN/m^2")
print(f"Undrained Shear Strength, c_uk = {c_uk:.2f} kN/m^2")
print(f"Shape Factor, s_c = {s_c:.2f}")
print(f"Inclination Factor, i_c = {i_c:.2f}")
print(f"Bearing Capacity Factor, (pi + 2) = {Nc:.4f}")
print(f"Partial Factor for Bearing Resistance, gamma_R = {gamma_R:.2f}\n")

print("Final Equation with numbers:")
print(f"R_d = ({A:.2f} * (({math.pi:.4f} + 2) * {c_uk:.2f} * {s_c:.2f} * {i_c:.2f} + {q:.2f})) / {gamma_R:.2f}")
print(f"R_d = ({A:.2f} * ({Nc:.4f} * {c_uk:.2f} * {s_c:.2f} * {i_c:.2f} + {q:.2f})) / {gamma_R:.2f}")
print(f"R_d = ({A:.2f} * ({q_gross_k:.2f})) / {gamma_R:.2f}")
print(f"R_d = {R_k:.2f} kN / {gamma_R:.2f}")

print("\n--- Final Answer ---")
print(f"The design resistance under combination 1 for undrained condition is: {R_d:.2f} kN")
<<<1349.32>>>