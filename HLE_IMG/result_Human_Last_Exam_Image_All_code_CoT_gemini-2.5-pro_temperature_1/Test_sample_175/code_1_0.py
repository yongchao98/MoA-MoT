import math

# --- Given Parameters ---
# Foundation geometry
B = 2.2  # Foundation width in m
L = 2.2  # Foundation length in m
D = 2.0  # Foundation depth in m

# Soil properties
c_uk = 57.0  # Characteristic undrained shear strength in kN/m^2 (kPa)

# Eurocode 7, DA1-1 Partial Factors
gamma_cu = 1.4  # Partial factor for undrained shear strength

# --- Calculation Steps ---

# 1. Effective Area (A')
A_prime = B * L

# 2. Shape Factor (s_c) for a square base
s_c = 1.2

# 3. Depth Factor (d_c)
d_c = 1 + 0.4 * (D / B)

# 4. Inclination Factor (i_c) for vertical load
i_c = 1.0

# 5. Design Undrained Shear Strength (c_ud)
c_ud = c_uk / gamma_cu

# 6. Design Net Bearing Capacity (q_net,d)
# q_net,d = (pi + 2) * c_ud * s_c * d_c * i_c
q_net_d = (math.pi + 2) * c_ud * s_c * d_c * i_c

# 7. Design Bearing Resistance (R_d)
R_d = q_net_d * A_prime

# --- Output the results ---
print("Calculation of Design Bearing Resistance (R_d) for Undrained Condition (EC7, DA1-1)")
print("---------------------------------------------------------------------------------")
print("The formula for design resistance is: R_d = A' * q_net,d")
print("where q_net,d = (π + 2) * c_ud * s_c * d_c * i_c")
print("and c_ud = c_uk / γ_cu\n")

print("Step 1: Calculate intermediate factors")
print(f"  - Foundation Area (A'): {B:.1f} m * {L:.1f} m = {A_prime:.2f} m^2")
print(f"  - Shape Factor (s_c): {s_c:.3f} (for square base)")
print(f"  - Depth Factor (d_c): 1 + 0.4 * ({D:.1f} / {B:.1f}) = {d_c:.3f}")
print(f"  - Inclination Factor (i_c): {i_c:.3f} (for vertical load)\n")

print("Step 2: Calculate design shear strength (c_ud)")
print(f"  - c_ud = c_uk / γ_cu = {c_uk:.1f} kN/m^2 / {gamma_cu:.1f} = {c_ud:.3f} kN/m^2\n")

print("Step 3: Calculate the design resistance (R_d)")
print("R_d = A' * (π + 2) * c_ud * s_c * d_c * i_c")
print(f"R_d = {A_prime:.2f} * ({math.pi:.3f} + 2) * {c_ud:.3f} * {s_c:.3f} * {d_c:.3f} * {i_c:.3f}")
print(f"R_d = {A_prime:.2f} * {q_net_d:.2f}")
print(f"R_d = {R_d:.2f} kN\n")

print("The design resistance under combination 1 for undrained condition is approximately {:.0f} kN.".format(R_d))

# Final answer format
print(f"\n<<<{R_d:.1f}>>>")