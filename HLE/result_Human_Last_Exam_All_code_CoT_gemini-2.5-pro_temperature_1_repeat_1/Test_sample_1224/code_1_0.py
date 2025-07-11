import math

# Titan 4-bit architecture uses fractions with integers from 0-15.
# We simulate this by choosing approximations that respect this constraint.

# --- Step 1: Define Constants and Approximations ---

# Drop height h (m)
h_val = 5000
h_num, h_den, h_exp = 5, 1, 3  # Represents 5/1 * 10^3

# Gravitational Constant G (m^3 kg^-1 s^-2)
# True value is ~6.674e-11. We approximate G as 2/3 * 1e-10 = 6.66...e-11
G_num, G_den, G_exp = 2, 3, -10

# Pi
# True value is ~3.14159. We approximate pi as 3/1.
pi_num, pi_den = 3, 1

# Pandora's shell density (kg/m^3)
# Value is 300. We represent as 3/1 * 10^2
rho_s_num, rho_s_den, rho_s_exp = 3, 1, 2

# Pandora's equatorial radius (m)
# Value is 2,000,000 m. A factor of 2 creates numbers > 15 in multiplication.
# We approximate it as 1,875,000 m = 15/8 * 1e6 to stay within limits.
r_p_num, r_p_den, r_p_exp = 15, 8, 6

# --- Step 2: Calculate g on Titan ---
# g = (4/3) * pi * G * rho_shell * r_planet
# We group multiplications to keep intermediate values small.
# (4/3 * pi) = 4/3 * 3/1 = 4/1
term1_num, term1_den = 4, 1

# (G * rho_shell) = (2/3e-10 * 3/1e2) = 2/1e-8
term2_num, term2_den, term2_exp = 2, 1, -8

# (term1 * term2) = 4/1 * 2/1e-8 = 8/1e-8
term3_num, term3_den, term3_exp = 8, 1, -8

# (term3 * r_planet) = 8/1e-8 * 15/8e6 = 15/1e-2
g_num, g_den, g_exp = 15, 1, -2
g_val = (g_num / g_den) * (10**g_exp)

# --- Step 3: Calculate 2h/g ---
# S = (2 * h) / g
h2_val = 2 * h_val
h2_exp = 4 # 1e4
h2_num, h2_den = 1,1

# S = (1/1 * 1e4) / (15/1 * 1e-2) = 1/15 * 1e6
S_num, S_den, S_exp = 1, 15, 6
S_val = (S_num / S_den) * (10**S_exp)

# --- Step 4: Calculate t = sqrt(S) ---
# We must use a pre-computed approximation for sqrt.
# t = sqrt(1/15 * 1e6) = sqrt(1/15) * 1e3
# The value sqrt(1/15) is ~0.258. The closest 4-bit fraction is 1/4.
sqrt_approx_num, sqrt_approx_den = 1, 4

# t = (1/4) * 1e3 = 250
t_calc_num = sqrt_approx_num
t_calc_den = sqrt_approx_den
t_calc_exp = 3
t_calculated = (t_calc_num / t_calc_den) * (10**t_calc_exp)

# --- Step 5: Final Equation Output ---
print("Final Calculation on Titan:")
print(f"h = {h_val}")
print(f"g ≈ {g_num}/{g_den}e{g_exp} = {g_val}")
final_equation = f"t = sqrt(2 * {h_val} / {g_val}) ≈ {t_calculated}"
print(final_equation)
print("\nEach number in the final equation:")
print(f"2 = {2}")
print(f"h = {h_val}")
print(f"g = {g_val}")
print(f"t = {t_calculated}")


# --- For Error Calculation (done on a standard computer) ---
# More precise values
G_true = 6.674e-11
pi_true = 3.1415926535
rho_s_true = 300
rho_c_true = 1200
r_p_true = 2e6
r_c_true = 1e5

# g_true = 4/3 * pi * G * [rho_s*r_p + (rho_c-rho_s)*r_c^3/r_p^2]
g_true_val = (4/3 * pi_true * G_true * 
             (rho_s_true * r_p_true + 
             (rho_c_true - rho_s_true) * r_c_true**3 / r_p_true**2))
t_true = math.sqrt(2 * h_val / g_true_val)
abs_error = abs(t_calculated - t_true)

# print(f"\n--- Error Analysis ---")
# print(f"True g: {g_true_val:.4f} m/s^2")
# print(f"True t: {t_true:.2f} s")
# print(f"Calculated t: {t_calculated:.2f} s")
# print(f"Smallest absolute error e: {abs_error:.2f}")

# The final answer format required by the user
# print(f"\n<<<Y[{abs_error:.2f}]>>>")
final_answer = f"Y[{abs_error:.2f}]"