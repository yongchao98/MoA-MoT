import math

# --- Given Parameters ---
# Soil properties from the problem description
c_uk = 57  # Characteristic undrained shear strength (kN/m^2)

# Foundation geometry
B = 2.2  # Foundation width (m)
L = 2.2  # Foundation length (m)
D = 2.0  # Foundation depth (m)

# Partial factor of safety for bearing resistance (EC7, DA1-1, undrained)
gamma_R = 1.4

# --- Calculations ---
# 1. Effective foundation area (A')
A_prime = B * L

# 2. Shape factor (s_c) for a square foundation
s_c = 1 + 0.2 * (B / L)

# 3. Depth factor (d_c)
d_c = 1 + 0.4 * (D / B)

# 4. Inclination factor (i_c) for vertical load
i_c = 1.0

# 5. Bearing capacity factor (N_c) for undrained analysis
Nc_factor = math.pi + 2

# 6. Characteristic net bearing resistance (R_k)
R_k = A_prime * (Nc_factor * c_uk * s_c * d_c * i_c)

# 7. Design bearing resistance (R_d)
R_d = R_k / gamma_R

# --- Output ---
print("Calculation of the Design Bearing Resistance (R_d) for Undrained Condition")
print("------------------------------------------------------------------------")
print("The formula for the design resistance is R_d = R_k / gamma_R, where:")
print("R_k = A' * ( (pi + 2) * c_uk * s_c * d_c * i_c )")
print("\nStep 1: Calculate the necessary factors:")
print(f"Foundation Area (A'): {B:.1f} m * {L:.1f} m = {A_prime:.2f} m^2")
print(f"Shape Factor (s_c): 1 + 0.2 * ({B:.1f}/{L:.1f}) = {s_c:.3f}")
print(f"Depth Factor (d_c): 1 + 0.4 * ({D:.1f}/{B:.1f}) = {d_c:.3f}")
print(f"Inclination Factor (i_c): {i_c:.1f} (for vertical load)")
print(f"Undrained Shear Strength (c_uk): {c_uk} kN/m^2")
print(f"Partial Factor for Bearing (gamma_R): {gamma_R:.1f}")

print("\nStep 2: Substitute these values into the equation for R_d:")
print("R_d = ( A' * ( (pi + 2) * c_uk * s_c * d_c * i_c ) ) / gamma_R")

# Print the final equation with all the numbers substituted
print(f"R_d = ( {A_prime:.2f} * ( ({math.pi:.3f} + 2) * {c_uk} * {s_c:.3f} * {d_c:.3f} * {i_c:.1f} ) ) / {gamma_R:.1f}")

# Print the intermediate and final results
print(f"\nCalculating the result:")
print(f"R_k = {R_k:.1f} kN")
print(f"R_d = {R_k:.1f} kN / {gamma_R:.1f}")
print(f"R_d = {R_d:.1f} kN")

print("\nThe final design resistance is {:.1f} kN.".format(R_d))

# Final answer in the required format
final_answer = round(R_d, 1)
print(f'<<<{final_answer}>>>')