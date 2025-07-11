# All numbers are represented as a tuple: (numerator, denominator, exponent)
# For example, 50kg is (5, 1, 1), representing (5/1) * 10^1

# --- Initial Values ---
m_probe = (5, 1, 1)    # 50 kg
v_initial = (3, 1, 2)  # 300 m/s
h_altitude = (5, 1, 3) # 5000 m

# --- Calculation of net deceleration force F_net = m * a ---

# 1. Calculate required acceleration 'a' from kinematics: a = v_i^2 / (2*h)
v_i_squared = (v_initial[0]**2, v_initial[1]**2, v_initial[2]*2) # (9, 1, 4)
two_h = (2 * h_altitude[0], h_altitude[1], h_altitude[2])       # (10, 1, 3)
# Normalize two_h to have a single digit numerator for scientific notation
two_h_norm = (1, 1, 4)
# a = v_i^2 / (2h)
a_frac_num = v_i_squared[0] * two_h_norm[1]
a_frac_den = v_i_squared[1] * two_h_norm[0]
a_exp = v_i_squared[2] - two_h_norm[2]
a = (a_frac_num, a_frac_den, a_exp) # (9, 1, 0) -> Represents 9 m/s^2

# 2. Calculate F_net = m * a
# F_net_frac = (m_probe[0]/m_probe[1]) * (a[0]/a[1]) = (5/1) * (9/1)
# The resulting numerator 5*9=45 exceeds the 5-bit limit of 31.
# ACTION: Sacrifice precision. Approximate a=9 with a_approx=6 to proceed.
a_approx = (6, 1, 0)

# Perform the valid multiplication
f_net_num = m_probe[0] * a_approx[0] # 5 * 6 = 30 (OK)
f_net_den = m_probe[1] * a_approx[1] # 1 * 1 = 1 (OK)
f_net_exp = m_probe[2] + a_approx[2]
f_net = (f_net_num, f_net_den, f_net_exp) # (30, 1, 1) -> Represents 300 N

# --- Calculation of gravitational force F_g = m * g ---

# 1. Calculate g ≈ (5/3)*10^-1 based on simplified planetary physics
g = (5, 3, -1)

# 2. Calculate F_g = m * g
# F_g_frac = (m_probe[0]/m_probe[1]) * (g[0]/g[1]) = (5/1) * (5/3)
f_g_num = m_probe[0] * g[0] # 5 * 5 = 25 (OK)
f_g_den = m_probe[1] * g[1] # 1 * 3 = 3 (OK)
f_g_exp = m_probe[2] + g[2]
f_g = (f_g_num, f_g_den, f_g_exp) # (25, 3, 0) -> Represents 8.333... N

# --- Final Calculation F_rocket = F_net + F_g ---

# To add, exponents must be aligned. Align to 10^1.
# F_net is (30, 1, 1)
# F_g is (25, 3, 0) -> 8.333... = 0.8333... * 10^1. Fraction for 0.8333... is 5/6.
f_g_aligned = (5, 6, 1)

# Add fractional parts: (30/1) + (5/6)
# Result is 180/6 + 5/6 = 185/6. Numerator 185 exceeds 31.
# ACTION: Sacrifice precision. Approximate the smaller term f_g_aligned.
# (5/6) is approximated as (1/1)
f_g_aligned_approx = (1, 1, 1)

# Perform the valid addition
# (30/1) + (1/1) = 31/1
f_rocket_num = (f_net[0] * f_g_aligned_approx[1]) + (f_g_aligned_approx[0] * f_net[1])
f_rocket_den = f_net[1] * f_g_aligned_approx[1]
f_rocket_exp = f_net[2]
f_rocket = (f_rocket_num, f_rocket_den, f_rocket_exp) # (31, 1, 1) -> Represents 310 N

# --- Print Final Equation ---
print("The final calculation for the landing force is:")
# F_net part
f_net_eq = f"(( {m_probe[0]} / {m_probe[1]} )*10^{m_probe[2]} * ( {a_approx[0]} / {a_approx[1]} ))"
# F_g part
f_g_eq = f"(( {m_probe[0]} / {m_probe[1]} )*10^{m_probe[2]} * ( {g[0]} / {g[1]} )*10^{g[2]})"
# Simplified terms
f_net_val = f"( {f_net[0]} / {f_net[1]} )*10^{f_net[2]}"
f_g_val = f"( {f_g[0]} / {f_g[1]} )*10^{f_g[2]}"
# Final result
f_rocket_val = f"( {f_rocket[0]} / {f_rocket[1]} )*10^{f_rocket[2]}"

print(f"Force = {f_net_eq} + {f_g_eq}")
print(f"      = {f_net_val} + {f_g_val}")
print(f"      ≈ ( {f_net[0]} / {f_net[1]} )*10^{f_net[2]} + ( {f_g_aligned_approx[0]} / {f_g_aligned_approx[1]} )*10^{f_g_aligned_approx[2]}")
print(f"      = {f_rocket_val}")

final_force_value = (f_rocket[0]/f_rocket[1]) * (10**f_rocket[2])
print(f"Final Calculated Force: {final_force_value} N")

# --- Error Calculation ---
true_g = 0.1677
true_a = 9
true_m = 50
true_force = true_m * (true_a + true_g) # 458.385
abs_error = abs(final_force_value - true_force)

# The question asks for Y[e] or N0 as the final answer
# This format should be at the very end.
# print(f"<<<Y[{abs_error:.3f}]>>>")