import math

# Step 1: Define given parameters from the problem
B = 2.2  # m, foundation width
L = 2.2  # m, foundation length
D = 2.0  # m, foundation depth
c_uk = 57.0  # kPa (kN/m^2), characteristic undrained shear strength
depth_gwl = 0.6  # m, groundwater level depth
gamma_soil_1 = 18.0  # kN/m^3, unit weight of clay above GWL
gamma_soil_2 = 20.0  # kN/m^3, unit weight of clay below GWL

# Step 2: Define Eurocode 7, DA1-1 partial factors
gamma_cu = 1.0  # Partial factor for undrained shear strength
gamma_gamma = 1.0  # Partial factor for soil unit weight
gamma_R_v = 1.0  # Partial factor for bearing resistance

# Step 3: Calculate Foundation Area (A')
A_prime = B * L

# Step 4: Calculate total overburden pressure (q) at foundation base
q = (depth_gwl * gamma_soil_1) + ((D - depth_gwl) * gamma_soil_2)

# Step 5: Calculate bearing capacity factors
s_c = 1 + 0.2 * (B / L)
d_c = 1 + 0.2 * (D / B)
i_c = 1.0  # Inclination factor for vertical load

# Step 6: Determine Design Soil Properties
# Design undrained shear strength
c_ud = c_uk / gamma_cu
# Design surcharge pressure
q_d_surcharge = q / gamma_gamma
# Bearing capacity factor for undrained condition (phi=0)
Nc = math.pi + 2

# Step 7: Calculate design gross bearing capacity (q_d)
q_d = Nc * c_ud * s_c * d_c * i_c + q_d_surcharge

# Step 8: Calculate design bearing resistance (R_d)
R_d = (A_prime * q_d) / gamma_R_v

# Print the detailed solution
print("Calculation of the Design Bearing Resistance (R_d)")
print("==================================================")
print("The analysis is for undrained conditions using Eurocode 7, Design Approach 1, Combination 1.\n")
print("The general formula for design bearing resistance is: R_d = A' * q_d / γ_R,v")
print("where q_d = (π + 2) * c_ud * s_c * d_c * i_c + q_d_surcharge\n")

print(f"1. Foundation Area (A'):")
print(f"   A' = B * L = {B:.1f} m * {L:.1f} m = {A_prime:.2f} m²\n")

print(f"2. Total Overburden Pressure at foundation base (q):")
print(f"   q = ({depth_gwl:.1f} m * {gamma_soil_1:.1f} kN/m³) + (({D:.1f} m - {depth_gwl:.1f} m) * {gamma_soil_2:.1f} kN/m³) = {q:.2f} kPa\n")

print(f"3. Design Soil Parameters (for DA1-1, γ_M=1.0):")
print(f"   Design undrained shear strength, c_ud = c_uk / γ_cu = {c_uk:.1f} / {gamma_cu:.1f} = {c_ud:.1f} kPa")
print(f"   Design surcharge pressure, q_d_surcharge = q / γ_γ = {q:.2f} / {gamma_gamma:.1f} = {q_d_surcharge:.2f} kPa\n")

print(f"4. Bearing Capacity Factors:")
print(f"   Shape factor, s_c = 1 + 0.2 * (B/L) = 1 + 0.2 * ({B:.1f}/{L:.1f}) = {s_c:.3f}")
print(f"   Depth factor, d_c = 1 + 0.2 * (D/B) = 1 + 0.2 * ({D:.1f}/{B:.1f}) = {d_c:.3f}")
print(f"   Inclination factor, i_c = {i_c:.3f} (for vertical load)")
print(f"   Bearing capacity factor, N_c = π + 2 = {Nc:.3f}\n")

print("5. Final Calculation of Design Resistance (R_d):")
print("The final equation with the calculated values is:")
print(f"R_d = {A_prime:.2f} * (({Nc:.3f}) * {c_ud:.1f} * {s_c:.3f} * {d_c:.3f} + {q_d_surcharge:.2f}) / {gamma_R_v:.1f}")

# Breaking down the calculation for clarity
q_d_net_part = Nc * c_ud * s_c * d_c * i_c
print(f"R_d = {A_prime:.2f} * ({q_d_net_part:.2f} + {q_d_surcharge:.2f}) / {gamma_R_v:.1f}")
print(f"R_d = {A_prime:.2f} * {q_d:.2f} / {gamma_R_v:.1f}\n")
print(f"The design resistance under combination 1 is R_d = {R_d:.1f} kN")
print("<<<2200.0>>>")