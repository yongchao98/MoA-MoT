import math

# Step 1: Define given parameters from the problem description.
B = 2.2  # Foundation width in m
L = 2.2  # Foundation length in m
D_f = 2.0  # Foundation depth in m
D_w = 0.6  # Groundwater level depth in m
gamma1 = 18.0  # Unit weight of clay above GWL in kN/m^3
gamma2 = 20.0  # Unit weight of clay below GWL in kN/m^3
c_uk = 57.0  # Characteristic undrained shear strength in kN/m^2 (kPa)
# Using a standard partial factor for bearing resistance as per Eurocode 7 (DA2, R2 set).
gamma_Rv = 1.4 

print("Step-by-step calculation for the design resistance:")
print("-" * 50)

# Step 2: Calculate the effective area of the foundation (A').
A_prime = B * L
print(f"1. Foundation Effective Area (A'):")
print(f"   A' = B * L = {B} m * {L} m = {A_prime:.2f} m^2\n")

# Step 3: Calculate the total overburden pressure (q) at the foundation base.
q = (D_w * gamma1) + ((D_f - D_w) * gamma2)
print(f"2. Overburden Pressure at Foundation Base (q):")
print(f"   q = (depth_above_gw * gamma_above) + (depth_below_gw * gamma_below)")
print(f"   q = ({D_w} m * {gamma1} kN/m^3) + (({D_f} m - {D_w} m) * {gamma2} kN/m^3)")
print(f"   q = {D_w * gamma1:.2f} kN/m^2 + {D_f - D_w:.2f} m * {gamma2} kN/m^3")
print(f"   q = {D_w * gamma1:.2f} kN/m^2 + {(D_f - D_w) * gamma2:.2f} kN/m^2 = {q:.2f} kN/m^2\n")

# Step 4: Calculate bearing capacity factors.
# Shape factor (s_c) for a square footing in undrained conditions.
s_c = 1 + 0.2 * (B / L)
print(f"3. Shape Factor (s_c):")
print(f"   For a square footing, s_c = 1 + 0.2 * (B/L) = 1 + 0.2 * ({B}/{L}) = {s_c:.2f}\n")

# Inclination factor (i_c) for vertical load.
i_c = 1.0
print(f"4. Inclination Factor (i_c):")
print(f"   For a vertical load, i_c = {i_c:.2f}\n")

# Step 5: Calculate characteristic bearing resistance (R_k).
# Formula: R_k = A' * [ (pi + 2) * c_uk * s_c * i_c + q ]
pi = math.pi
net_bearing_pressure = (pi + 2) * c_uk * s_c * i_c
total_bearing_pressure = net_bearing_pressure + q
R_k = A_prime * total_bearing_pressure

print(f"5. Characteristic Bearing Resistance (R_k):")
print(f"   R_k = A' * [ (pi + 2) * c_uk * s_c * i_c + q ]")
print(f"   R_k = {A_prime:.2f} m^2 * [ ({pi:.3f} + 2) * {c_uk} kN/m^2 * {s_c:.2f} * {i_c:.2f} + {q:.2f} kN/m^2 ]")
print(f"   R_k = {A_prime:.2f} * [ {net_bearing_pressure:.2f} kN/m^2 + {q:.2f} kN/m^2 ]")
print(f"   R_k = {A_prime:.2f} * [ {total_bearing_pressure:.2f} kN/m^2 ] = {R_k:.2f} kN\n")

# Step 6: Calculate the design resistance (R_d).
R_d = R_k / gamma_Rv
print(f"6. Final Design Resistance (R_d):")
print(f"   R_d = R_k / γ_Rv")
print(f"   R_d = {R_k:.2f} kN / {gamma_Rv}")
print(f"   R_d = {R_d:.2f} kN")
print("-" * 50)
print(f"\nThe design resistance under combination 1 for undrained condition is {R_d:.1f} kN.")
