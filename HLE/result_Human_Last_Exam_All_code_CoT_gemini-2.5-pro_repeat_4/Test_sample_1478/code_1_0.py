# Titan 6-bit Computer Program for Gravitational Force Calculation

print("--- Titan 6-bit Architecture: Pandora Black Hole Gravity Calculation ---")
print("This program simulates the calculation of the gravitational force from the exoplanet Pandora, were it a black hole, on a 50kg probe.")
print("All calculations adhere to the 6-bit (0-63) fractional arithmetic constraints of the Titan architecture.")

# --- Step 0: Define Constants with 6-bit Fractional Approximations ---
# Gravitational Constant, G ≈ 6.674e-11 N⋅m²/kg². Approximated as 20/3 ≈ 6.667.
G_n, G_d, G_e = 20, 3, -11
# Pi, π ≈ 3.14159. Approximated as 22/7 ≈ 3.1428.
pi_n, pi_d = 22, 7
# Pandora Radius, r = 2000 km = 2e6 m.
r_n, r_d, r_e = 2, 1, 6
# Pandora Density, ρ = 1.2 tons/m³ = 1200 kg/m³.
rho_n, rho_d, rho_e = 12, 1, 2
# Probe Mass, m = 50 kg.
m_n, m_d = 50, 1
# Speed of Light, c ≈ 3e8 m/s.
c_n, c_d, c_e = 3, 1, 8
# Distance from event horizon, d_h = 1 km = 1000 m.
dh_n, dh_d, dh_e = 1, 1, 3

print("\n--- Step 1: Calculate Mass (M) of Pandora ---")
print(f"M = 4/3 * π * r³ * ρ")
# M = (4/3) * (22/7) * (2e6 m)³ * (1200 kg/m³)
# First, calculate the coefficient for the volume of a sphere: 4/3 * π
c1_n = 4 * pi_n  # 4 * 22 = 88
c1_d = 3 * pi_d  # 3 * 7 = 21
print(f"Intermediate calculation for 4/3 * π: {c1_n}/{c1_d}.")
print(f"Constraint Violation: Numerator {c1_n} is > 63.")
# Simplify 88/21 ≈ 4.19. A good 6-bit approximation is 25/6 ≈ 4.167
V_coeff_n, V_coeff_d = 25, 6
print(f"Simplifying to {V_coeff_n}/{V_coeff_d} to maintain 6-bit constraint.\n")

# Now calculate the full mass mantissa
# M_mantissa = (25/6) * (r_n)³ * rho_n = 25/6 * 2³ * 12
# We can simplify before the multiplication would overflow: (12/6) = 2
M_mant_n_simplified = V_coeff_n * (r_n**3) * (rho_n // V_coeff_d) # 25 * 8 * 2 = 400
print(f"Intermediate mass mantissa calculation: {V_coeff_n}/{V_coeff_d} * {r_n**3} * {rho_n} = {M_mant_n_simplified}/1.")
print(f"Constraint Violation: Numerator {M_mant_n_simplified} is > 63.")
# Simplify 400 by expressing it in scientific notation: 4 * 10^2
M_n, M_d = 4, 1
exp_adj = 2
M_e = (r_e * 3) + rho_e + exp_adj # 6*3 + 2 + 2 = 22
print(f"Simplifying to {M_n}/{M_d} and adjusting exponent. M ≈ {M_n}/{M_d} x 10^{M_e} kg.\n")


print("--- Step 2: Calculate Schwarzschild Radius (Rs) ---")
print("Rs = 2 * G * M / c²")
# c² = (3e8)² = 9e16
c2_n, c2_d, c2_e = c_n**2, c_d**2, c_e*2

# Rs numerator: 2 * G * M
Rs_num_n = 2 * G_n * M_n # 2 * 20 * 4 = 160
Rs_num_d = G_d * M_d # 3 * 1 = 3
Rs_num_e = G_e + M_e # -11 + 22 = 11
print(f"Intermediate Rs numerator mantissa: {Rs_num_n}/{Rs_num_d}.")
print(f"Constraint Violation: Numerator {Rs_num_n} is > 63.")
# Simplify 160/3 ≈ 53.33. A good 6-bit approximation is 53/1
Rs_num_n_s, Rs_num_d_s = 53, 1
print(f"Simplifying to {Rs_num_n_s}/{Rs_num_d_s} to maintain 6-bit constraint.\n")

# Final Rs calculation by dividing by c²
Rs_n = Rs_num_n_s * c2_d # 53 * 1
Rs_d = Rs_num_d_s * c2_n # 1 * 9
Rs_e = Rs_num_e - c2_e # 11 - 16 = -5
print(f"Calculated Schwarzschild Radius: Rs ≈ {Rs_n}/{Rs_d} x 10^{Rs_e} m.\n")

print("--- Step 3: Calculate Total Distance (d) from Center ---")
print("d = Rs + distance_from_horizon (1000 m)")
# d = 53/9e-5 m + 1000 m
# The value of Rs (≈ 5.9e-5 m) is negligible compared to 1000 m.
print(f"The term Rs ({Rs_n}/{Rs_d}e{Rs_e}) is negligible compared to the probe's distance (1000 m).")
d_n, d_d, d_e = dh_n, dh_d, dh_e
d2_n, d2_d, d2_e = d_n**2, d_d**2, d_e*2
print(f"Simplification: d is approximated as {d_n}/{d_d} x 10^{d_e} m. Thus, d² = {d2_n}/{d2_d} x 10^{d2_e}.\n")

print("--- Step 4: Calculate Final Gravitational Force (F) ---")
print("F = G * M * m / d²")
# Numerator: G * M * m
F_num_n = G_n * M_n * m_n # 20 * 4 * 50 = 4000
F_num_d = G_d * M_d * m_d # 3 * 1 * 1 = 3
F_num_e = G_e + M_e # -11 + 22 = 11
print(f"Intermediate force numerator mantissa: {F_num_n}/{F_num_d}.")
print(f"Constraint Violation: Numerator {F_num_n} is > 63.")
# Simplify 4000/3 by expressing 4000 as 4 * 10^3. The mantissa becomes 4/3.
F_num_n_s, F_num_d_s = 4, 3
exp_adj_f = 3
F_num_e_s = F_num_e + exp_adj_f # 11 + 3 = 14
print(f"Simplifying to {F_num_n_s}/{F_num_d_s} and adjusting exponent.\n")

# Final Force calculation
F_n = F_num_n_s * d2_d # 4 * 1 = 4
F_d = F_num_d_s * d2_n # 3 * 1 = 3
F_e = F_num_e_s - d2_e # 14 - 6 = 8

print("--- Final Result ---")
print("The final equation with all approximated values is:")
print(f"F = (({G_n}/{G_d} x 10^({G_e})) * ({M_n}/{M_d} x 10^({M_e})) * ({m_n}/{m_d})) / ({d_n}/{d_d} x 10^({d_e}))²")
print(f"\nFinal Calculated Force: F = {F_n}/{F_d} x 10^{F_e} N")