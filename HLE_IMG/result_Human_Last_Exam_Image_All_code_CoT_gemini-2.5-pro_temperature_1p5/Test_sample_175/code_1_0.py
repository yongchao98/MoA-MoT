import math

# Step 1: Define the given parameters
B = 2.2  # Foundation width in m
L = 2.2  # Foundation length in m (B=L)
d = 2.0  # Foundation depth in m
d_w = 0.6  # Groundwater level depth in m
c_uk = 57.0  # Characteristic undrained shear strength in kN/m^2
gamma_soil_above = 18.0  # Unit weight of clay above water table in kN/m^3
gamma_soil_below = 20.0  # Unit weight of clay below water table in kN/m^3
gamma_R_v = 1.4  # Partial factor for bearing resistance for DA1 Combination 1

print("Step 1: Calculate the effective area of the foundation (A').")
A_prime = B * L
print(f"A' = B * L = {B} m * {L} m = {A_prime:.2f} m^2\n")

print("Step 2: Calculate the total overburden pressure (q) at the foundation base.")
q = (d_w * gamma_soil_above) + ((d - d_w) * gamma_soil_below)
print(f"q = (d_w * γ_soil_above) + ((d - d_w) * γ_soil_below)")
print(f"q = ({d_w} m * {gamma_soil_above} kN/m^3) + (({d} m - {d_w} m) * {gamma_soil_below} kN/m^3)")
print(f"q = {d_w * gamma_soil_above:.2f} kN/m^2 + {(d - d_w):.1f} m * {gamma_soil_below} kN/m^3 = {d_w * gamma_soil_above:.2f} + {(d - d_w) * gamma_soil_below:.2f} = {q:.2f} kN/m^2\n")

print("Step 3: Determine the bearing capacity, shape, inclination and base factors.")
# For undrained conditions (phi_u = 0):
N_c = math.pi + 2
# For a square footing (B/L = 1):
s_c = 1.2
# For a vertical load:
i_c = 1.0
# For a horizontal base:
b_c = 1.0
print(f"Bearing capacity factor (undrained): N_c = π + 2 = {N_c:.3f}")
print(f"Shape factor (square): s_c = {s_c}")
print(f"Load inclination factor (vertical load): i_c = {i_c}")
print(f"Base inclination factor (horizontal base): b_c = {b_c}\n")

print("Step 4: Calculate the characteristic bearing resistance (R_k).")
R_k = A_prime * (N_c * c_uk * s_c * i_c * b_c + q)
print(f"R_k = A' * (N_c * c_uk * s_c * i_c * b_c + q)")
print(f"R_k = {A_prime:.2f} * (({N_c:.3f} * {c_uk} * {s_c} * {i_c} * {b_c}) + {q:.2f})")
term1 = N_c * c_uk * s_c * i_c * b_c
print(f"R_k = {A_prime:.2f} * ({term1:.2f} + {q:.2f})")
print(f"R_k = {A_prime:.2f} * {term1 + q:.2f} = {R_k:.2f} kN\n")

print("Step 5: Calculate the design resistance (R_d) for undrained conditions.")
R_d = R_k / gamma_R_v
print(f"The partial factor of safety for bearing resistance, γ_R,v = {gamma_R_v}")
print(f"R_d = R_k / γ_R,v")
print(f"R_d = {R_k:.2f} kN / {gamma_R_v} = {R_d:.2f} kN")

print("\nThe final answer is the design resistance R_d.")
print(f"Final Answer = {R_d:.2f} kN")
<<<1350.98>>>