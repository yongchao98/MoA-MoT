import math

# Helper class to represent a number in Titan's architecture
# A number is represented as (numerator/denominator) * 10^exponent
class TitanNumber:
    def __init__(self, num, den, exp=0):
        if not (0 <= num <= 31 and 1 <= den <= 31 and -16 <= exp <= 15):
            # In a real Titan system, this would raise an error.
            # For this simulation, we assume valid inputs are chosen.
            pass
        self.num = num
        self.den = den
        self.exp = exp

    def __str__(self):
        if self.exp == 0:
            return f"({self.num}/{self.den})"
        return f"({self.num}/{self.den}) * 10^{self.exp}"

    def to_float(self):
        return (self.num / self.den) * (10 ** self.exp)

# --- Step 1: Define constants and perform real-world calculation for target value ---
m_probe_real = 50.0  # kg
v_i_real = 300.0  # m/s
d_real = 5000.0  # m
G_real = 6.674e-11 # N m^2 / kg^2
# Mass of Pandora (calculated from volumes and densities)
# V_core = 4/3*pi*(1e5)^3, V_shell = 4/3*pi*((2e6)^3 - (1e5)^3)
# M_core = V_core * 1200, M_shell = V_shell * 300
M_pandora_real = 1.0057e22 # kg
# Distance from center when landing system activates (equatorial radius + altitude)
r_real = 2e6 + 5000 # m

# Real-world force calculation
a_d_real = v_i_real**2 / (2 * d_real)
f_d_real = m_probe_real * a_d_real
f_g_real = (G_real * M_pandora_real * m_probe_real) / (r_real**2)
f_total_real = f_d_real + f_g_real

# --- Step 2: Perform calculation using Titan architecture ---

print("--- Titan Calculation Steps ---")

# Represent initial values as Titan numbers
# To calculate F_d = m * a = m * (v_i^2 / (2d)), we need to choose representations
# that allow multiplication without overflow (num1*num2 <= 31).
# m = 50 = (1/2)*10^2. This makes its numerator 1, which is easy to multiply.
m_probe = TitanNumber(1, 2, exp=2)
v_i = TitanNumber(3, 1, exp=2)
d = TitanNumber(5, 1, exp=3)

print(f"Probe Mass (m) = {m_probe_real} kg -> {m_probe}")
print(f"Initial Velocity (v_i) = {v_i_real} m/s -> {v_i}")
print(f"Distance (d) = {d_real} m -> {d}")
print("-" * 20)

# Calculate deceleration force, F_d
# a_d = v_i^2 / (2*d)
v_i_squared = TitanNumber(v_i.num**2, v_i.den**2, exp=v_i.exp*2) # (9/1)*10^4
two_d = TitanNumber(2 * d.num, d.den, exp=d.exp) # (10/1)*10^3 = (1/1)*10^4
# To divide, we multiply by the reciprocal and subtract exponents.
# (9/1) / (10/1) -> (9/1) * (1/10). Numerators are too large.
# Let's simplify 2d first: (10/1)*10^3 -> (1/1)*10^4
two_d_simplified = TitanNumber(1, 1, exp=4)
a_d = TitanNumber(v_i_squared.num * two_d_simplified.den, v_i_squared.den * two_d_simplified.num, exp=v_i_squared.exp - two_d_simplified.exp)
# a_d = (9*1 / 1*1) * 10^(4-4) = (9/1)
print(f"Calculated acceleration (a_d) = {a_d}")

# F_d = m * a_d
# F_d = (1/2)*10^2 * (9/1) = (1*9 / 2*1) * 10^2 = (9/2)*10^2
f_d = TitanNumber(m_probe.num * a_d.num, m_probe.den * a_d.den, exp=m_probe.exp + a_d.exp)
print(f"Deceleration Force (F_d) = {f_d} = {f_d.to_float()} N")
print("-" * 20)

# Calculate gravitational force, F_g
# F_g = G * M * m / r^2
G = TitanNumber(20, 3, exp=-11) # Approx 6.67e-11
M = TitanNumber(1, 1, exp=22)   # Approx 1.0057e22
r_squared = TitanNumber(4, 1, exp=12) # Approx (2.005e6)^2 = 4.02e12

# F_g = (G * M * m) / r^2
# G*M = (20/3)*10^-11 * (1/1)*10^22 = (20/3)*10^11
gm = TitanNumber(G.num * M.num, G.den * M.den, exp=G.exp + M.exp)
# G*M*m = (20/3)*10^11 * (1/2)*10^2 = (20/6)*10^13 = (10/3)*10^13
gmm = TitanNumber(gm.num * m_probe.num, gm.den * m_probe.den, exp=gm.exp + m_probe.exp)
gmm_simplified = TitanNumber(10, 3, exp=13) # Simplify 20/6 -> 10/3
# F_g = (10/3)*10^13 / (4/1)*10^12 = (10/12)*10^1 = (5/6)*10^1
f_g = TitanNumber(gmm_simplified.num * r_squared.den, gmm_simplified.den * r_squared.num, exp=gmm_simplified.exp - r_squared.exp)
# Simplify 10/12 -> 5/6
f_g_simplified = TitanNumber(5, 6, exp=1)
print(f"Gravitational Force (F_g) = {f_g_simplified} = {f_g_simplified.to_float():.3f} N")
print("-" * 20)

# Add the forces: F_total = F_d + F_g
# F_total = (9/2)*10^2 + (5/6)*10^1 = 450 + 8.333...
# To add, exponents must match. Let's use exp=2.
# F_g = (5/6)*10^1 = (5/60)*10^2 = (1/12)*10^2
f_g_exp2 = TitanNumber(1, 12, exp=2)
# F_total = (9/2)*10^2 + (1/12)*10^2 = (9/2 + 1/12)*10^2
# 9/2 + 1/12 = (54/12 + 1/12) = 55/12.
# The numerator 55 is > 31, so this operation is illegal.
# We must approximate the result.
# The target value is 458.348... = 4.58348 * 10^2
# We need to find the best fraction num/den <= 31 for 4.58348
# Possible fractions: 9/2=4.5, 14/3=4.66, 23/5=4.6
# Error for 23/5: |4.6 - 4.58348| = 0.01652
# Error for 9/2: |4.5 - 4.58348| = 0.08348
# Error for 14/3: |4.667 - 4.58348| = 0.08352
# The best approximation is 23/5.
f_total_approx = TitanNumber(23, 5, exp=2)
print(f"Addition F_d + F_g results in overflow (55/12).")
print(f"Approximating the result to the nearest valid Titan number.")
print(f"Final Total Force (F_total) ≈ {f_total_approx} = {f_total_approx.to_float()} N")
print("-" * 20)

# --- Step 3: Final Equation and Error Calculation ---
error = abs(f_total_approx.to_float() - f_total_real)

print("Final Equation:")
# The problem asks to output each number in the final equation.
# The equation is F_total = F_d + F_g. We show the final approximated result
# and the calculated components.
print(f"{f_total_approx} = {f_d} + {f_g_simplified}")
print("\nNote: The above equation is an approximation, as the addition required sacrificing precision to meet the 5-bit constraint.")

print("\n--- Result ---")
print(f"Real-world force: {f_total_real:.3f} N")
print(f"Titan-computed force: {f_total_approx.to_float()} N")
print(f"Smallest absolute error: {error:.3f} N")

# Final answer format
final_answer = f"Y[{error:.3f}]"
print(f"\nFinal Answer in specified format: {final_answer}")

# The final answer is Y[1.652]
# The problem asks to return the answer directly in the specified format.
# The code above serves as the justification for the answer.
# <<<Y[1.652]>>>