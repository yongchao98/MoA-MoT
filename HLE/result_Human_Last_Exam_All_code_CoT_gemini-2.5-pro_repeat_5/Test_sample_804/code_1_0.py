import math

# Helper class for Titan fractions.
# In a real scenario, this would be hardware logic.
# For this simulation, it helps format the output.
class TitanFraction:
    def __init__(self, num, den=1, name=""):
        # We assume inputs are valid for this simulation
        self.num = num
        self.den = den
        self.name = name

    def __str__(self):
        return f"({self.num}/{self.den})"

# 1. Start with the simplified physics formula.
# The core's mass and probe's altitude are considered negligible
# to keep calculations within Titan's constraints.
# F ≈ G * (4/3)π * ρ_shell * m_probe * r_total
print("Step 1: Simplify the physics formula.")
print("F = G * (4/3) * π * ρ_shell * m_probe * r_total")
print("-" * 20)

# 2. Group terms and approximate them for Titan's architecture.
print("Step 2: Group terms and create fractional approximations for Titan.")

# Group 1: K = G * (4/3) * π
# Real value ≈ 2.796e-10. Approximate as 0.25e-9 = 1/4 * 10^-9.
K_approx = TitanFraction(1, 4, name="G * (4/3) * π")
print(f"Approximation for {K_approx.name}: {K_approx} * 10^-9")

# Group 2: Dm = ρ_shell * m_probe
# Real value = 300 * 30 = 9000. Approximate as 10000 = 1 * 10^4.
Dm_approx = TitanFraction(1, 1, name="ρ_shell * m_probe")
print(f"Approximation for {Dm_approx.name}: {Dm_approx} * 10^4")

# Group 3: Rt = r_total
# Real value = 1,000,000 = 1 * 10^6.
Rt_approx = TitanFraction(1, 1, name="r_total")
print(f"Approximation for {Rt_approx.name}: {Rt_approx} * 10^6")
print("-" * 20)

# 3. Perform the calculation using fractional arithmetic.
print("Step 3: Calculate the force using fractional arithmetic.")
# The coefficients are multiplied first.
f_coeff_num = K_approx.num * Dm_approx.num * Rt_approx.num
f_coeff_den = K_approx.den * Dm_approx.den * Rt_approx.den
f_coeff = TitanFraction(f_coeff_num, f_coeff_den)

# The exponents of 10 are summed up.
k_exp = -9
dm_exp = 4
rt_exp = 6
f_exp = k_exp + dm_exp + rt_exp

print(f"F_coeff = {K_approx} * {Dm_approx} * {Rt_approx} = {f_coeff}")
print(f"F_exp = ({k_exp}) + ({dm_exp}) + ({rt_exp}) = {f_exp}")
print("-" * 20)

# 4. Combine coefficient and exponent for the final result.
print("Step 4: Combine parts to get the final result.")
final_num = f_coeff.num * (10**f_exp)
final_den = f_coeff.den

# Simplify the final fraction
common_divisor = math.gcd(final_num, final_den)
final_frac = TitanFraction(final_num // common_divisor, final_den // common_divisor)

print(f"Final F = {f_coeff} * 10^{f_exp} = ({f_coeff.num}/{f_coeff.den}) * {10**f_exp} = {final_num}/{final_den}")
print(f"Simplified Final F = {final_frac.num} / {final_frac.den}")
print("-" * 20)

# 5. State the final answer and calculate the error.
real_F = 2.514  # More precise value
titan_F = final_frac.num / final_frac.den
error = abs(titan_F - real_F)

print(f"Can Titan be used? Yes.")
print(f"What is the smallest absolute error? e = |{titan_F} - {real_F}| = {error:.3f}")
print("Final Answer in required format:")
print(f"<<<Y[{error:.3f}]>>>")
