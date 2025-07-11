import math

# Step 1: Define given parameters from the problem
c_uk = 57.0  # Characteristic undrained shear strength (kN/m^2)
B = 2.2      # Foundation width (m)
L = 2.2      # Foundation length (m)
D_f = 2.0    # Foundation depth (m)
h_gw = 0.6   # Groundwater level depth (m)
gamma_dry = 18.0 # Unit weight of clay above groundwater (kN/m^3)
gamma_sat = 20.0 # Unit weight of clay below groundwater (kN/m^3)
gamma_R = 1.4  # Partial factor for bearing resistance (undrained)

# Step 2: Calculate foundation area and overburden pressure (q)
A_prime = B * L
q = (gamma_dry * h_gw) + (gamma_sat * (D_f - h_gw))

print("Step 1: Calculate Overburden Pressure (q) at foundation base")
print(f"q = (γ_dry * h_gw) + (γ_sat * (D_f - h_gw))")
print(f"q = ({gamma_dry} * {h_gw}) + ({gamma_sat} * ({D_f} - {h_gw}))")
print(f"q = {gamma_dry * h_gw:.2f} + {gamma_sat * (D_f - h_gw):.2f} = {q:.2f} kN/m^2\n")

# Step 3: Calculate bearing capacity factors for undrained condition (phi_u = 0)
Nc_term = math.pi + 2
s_c = 1.2  # Shape factor for square foundation
i_c = 1.0  # Inclination factor for vertical load
d_c = 1 + 0.4 * (D_f / B)

print("Step 2: Calculate Bearing Capacity Factors")
print(f"Shape factor for square foundation, s_c = {s_c}")
print(f"Depth factor, d_c = 1 + 0.4 * (D_f / B) = 1 + 0.4 * ({D_f} / {B}) = {d_c:.4f}")
print(f"Inclination factor for vertical load, i_c = {i_c}")
print(f"Bearing capacity factor, (π + 2) = {Nc_term:.4f}\n")

# Step 4: Calculate characteristic gross bearing capacity (q_uk)
q_uk = (Nc_term * c_uk * s_c * i_c * d_c) + q

print("Step 3: Calculate Characteristic Gross Bearing Capacity (q_uk)")
print(f"q_uk = ((π + 2) * c_uk * s_c * i_c * d_c) + q")
print(f"q_uk = (({Nc_term:.4f}) * {c_uk} * {s_c} * {i_c} * {d_c:.4f}) + {q:.2f}")
print(f"q_uk = {Nc_term * c_uk * s_c * i_c * d_c:.2f} + {q:.2f} = {q_uk:.2f} kN/m^2\n")

# Step 5: Calculate characteristic gross bearing resistance (R_k)
R_k = q_uk * A_prime

print("Step 4: Calculate Characteristic Gross Bearing Resistance (R_k)")
print(f"Foundation Area, A' = B * L = {B} * {L} = {A_prime:.2f} m^2")
print(f"R_k = q_uk * A'")
print(f"R_k = {q_uk:.2f} * {A_prime:.2f} = {R_k:.2f} kN\n")

# Step 6: Calculate the final design bearing resistance (R_d)
R_d = R_k / gamma_R

print("Step 5: Calculate Design Bearing Resistance (R_d)")
print(f"Partial factor for bearing resistance, γ_R = {gamma_R}")
print(f"R_d = R_k / γ_R")
print(f"R_d = {R_k:.2f} / {gamma_R} = {R_d:.2f} kN")

print(f"\nThe final design resistance is {R_d:.1f} kN.")

# Final answer in the specified format
print(f"\n<<<{R_d:.1f}>>>")