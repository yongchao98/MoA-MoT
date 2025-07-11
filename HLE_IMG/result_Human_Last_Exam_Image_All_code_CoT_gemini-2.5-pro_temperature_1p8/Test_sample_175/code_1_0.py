import math

# Step 1: Define given parameters from the problem statement
# Characteristic undrained shear strength (kN/m^2)
c_uk = 57.0
# Foundation width (m)
B = 2.2
# Foundation length (m)
L = 2.2
# Foundation depth (m)
D = 2.0
# For a square foundation under vertical load
i_c = 1.0 # Load inclination factor
b_c = 1.0 # Base inclination factor

# Step 2: Determine partial factors for Eurocode 7, Design Approach 1, Combination 1
gamma_cu = 1.0 # Partial factor for undrained shear strength
gamma_Rv = 1.0 # Partial factor for bearing resistance

# Step 3: Calculate foundation area and other factors
# Foundation Area A' (m^2)
A_prime = B * L
# Shape factor s_c for a square foundation in undrained analysis
s_c = 1.2
# Depth factor d_c for undrained analysis
d_c = 1 + 0.4 * (D / B)

# Step 4: Calculate design net bearing capacity (q_n,d)
# Design undrained shear strength (c_ud = c_uk / gamma_cu)
c_ud = c_uk / gamma_cu
# Bearing capacity factor for undrained conditions (Nc = pi + 2)
Nc = math.pi + 2
# Design net bearing capacity (kN/m^2)
q_n_d = Nc * c_ud * s_c * d_c * i_c * b_c

# Step 5: Calculate design bearing resistance (R_d)
# Design bearing resistance (kN)
R_d = (A_prime * q_n_d) / gamma_Rv

# Step 6: Output the results and the final calculation
print("--- Calculation of Design Bearing Resistance (Undrained, Combination 1) ---")
print(f"Foundation Width (B): {B} m")
print(f"Foundation Length (L): {L} m")
print(f"Foundation Area (A'): {A_prime:.2f} m^2")
print(f"Foundation Depth (D): {D} m")
print(f"Characteristic Undrained Shear Strength (c_uk): {c_uk} kN/m^2\n")

print("--- Calculation Factors ---")
print(f"Design Undrained Shear Strength (c_ud = c_uk / gamma_cu): {c_ud:.2f} kN/m^2")
print(f"Shape Factor (s_c): {s_c}")
print(f"Depth Factor (d_c = 1 + 0.4 * D/B): {d_c:.4f}")
print(f"Load Inclination Factor (i_c): {i_c}")
print(f"Base Inclination Factor (b_c): {b_c}")
print(f"Partial Resistance Factor (gamma_Rv): {gamma_Rv}\n")

print("--- Final Resistance Calculation ---")
print("The design resistance, Rd, is calculated using the formula:")
print("Rd = (A' * ( (pi + 2) * c_ud * s_c * d_c * i_c * b_c )) / gamma_Rv\n")

print("Substituting the values into the equation:")
print(f"Rd = ({A_prime:.2f} * (({math.pi:.4f} + 2) * {c_ud:.1f} * {s_c:.1f} * {d_c:.4f} * {i_c:.1f} * {b_c:.1f})) / {gamma_Rv:.1f}")
print(f"Rd = ({A_prime:.2f} * {q_n_d:.2f}) / {gamma_Rv:.1f}")
print(f"\nFinal Result:")
print(f"The design resistance under combination 1 for undrained condition is {R_d:.2f} kN.")
