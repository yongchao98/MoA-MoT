import math

# --- Given Parameters ---
B = 2.2  # Foundation width (m)
L = 2.2  # Foundation length (m)
d_foundation = 2.0  # Foundation depth (m)
d_water = 0.6  # Groundwater level depth (m)
c_uk = 57  # Undrained shear strength of clay (kN/m^2 or kPa)
gamma_clay_above = 18  # Unit weight of clay above groundwater (kN/m^3)
gamma_clay_below = 20  # Unit weight of clay below groundwater (kN/m^3)
gamma_Rv = 1.4  # Partial factor for bearing resistance (undrained, DA1)

# --- Calculation ---
print("Calculation of the Design Bearing Resistance (R_d)")
print("-------------------------------------------------")

# Step 1: Calculate foundation base area (A')
A_prime = B * L
print(f"1. Foundation base area (A'):")
print(f"   A' = B * L = {B} m * {L} m = {A_prime:.2f} m^2\n")

# Step 2: Calculate total overburden pressure (q) at foundation base
q = (gamma_clay_above * d_water) + (gamma_clay_below * (d_foundation - d_water))
print(f"2. Total overburden pressure at foundation base (q):")
print(f"   q = (γ_clay_above * d_water) + (γ_clay_below * (d_foundation - d_water))")
print(f"   q = ({gamma_clay_above} * {d_water}) + ({gamma_clay_below} * ({d_foundation} - {d_water}))")
print(f"   q = {gamma_clay_above * d_water:.2f} + {gamma_clay_below * (d_foundation - d_water):.2f} = {q:.2f} kPa\n")

# Step 3: Calculate characteristic gross bearing resistance (R_k)
# For a square footing, shape factor s_c = 1.2
# For a vertical load, inclination factor i_c = 1.0
pi = math.pi
s_c = 1.2
i_c = 1.0

# Calculate characteristic gross bearing pressure (q_gross_k)
q_gross_k = ((pi + 2) * c_uk * s_c * i_c) + q

print(f"3. Characteristic gross bearing resistance (R_k):")
print(f"   The formula is R_k = A' * ( (π + 2) * c_uk * s_c * i_c + q )")
print(f"   - Undrained bearing capacity factor (Nc): π + 2 = {pi + 2:.4f}")
print(f"   - Shape factor (s_c) for square footing = {s_c}")
print(f"   - Inclination factor (i_c) for vertical load = {i_c}")

# Calculate R_k
R_k = q_gross_k * A_prime

print(f"\n   First, we find the gross bearing pressure (q_gross_k):")
print(f"   q_gross_k = (({pi:.4f} + 2) * {c_uk} * {s_c} * {i_c}) + {q:.2f} = {q_gross_k:.2f} kPa")
print(f"\n   Then, R_k = q_gross_k * A'")
print(f"   R_k = {q_gross_k:.2f} kPa * {A_prime:.2f} m^2 = {R_k:.2f} kN\n")


# Step 4: Calculate the design resistance (R_d)
R_d = R_k / gamma_Rv

print(f"4. Design resistance (R_d):")
print(f"   The partial factor for undrained bearing resistance (γ_R,v) is {gamma_Rv}.")
print(f"   R_d = R_k / γ_R,v")
print(f"   R_d = {R_k:.2f} kN / {gamma_Rv} = {R_d:.2f} kN\n")

# Final Answer Summary
print("-------------------------------------------------")
print(f"The design resistance under combination 1 for undrained condition is {R_d:.2f} kN.")
print("-------------------------------------------------")
