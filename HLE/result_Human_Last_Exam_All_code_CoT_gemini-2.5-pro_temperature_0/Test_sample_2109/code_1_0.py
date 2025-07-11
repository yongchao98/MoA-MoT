from fractions import Fraction

# Step 1: Determine the optimal value of the parameter 'a'.
# From the analysis, the physical constraint Φ(z) >= 0 for z in [0, 3]
# leads to the condition a <= 1/32.
# The total energy E_total(a) is a decreasing function of 'a'.
# Therefore, the minimum energy occurs at the maximum possible value of 'a'.
a = Fraction(1, 32)

# Step 2: Determine coefficients 'b' and 'c' based on 'a'.
# From the fractional derivative constraints, we derived:
# 8a + 6b = 3  => b = (3 - 8a) / 6
# 32a + c = 1  => c = 1 - 32a
b = (3 - 8 * a) / 6
c = 1 - 32 * a

print(f"Derived coefficients for Φ(z) = az^2 + bz + c:")
print(f"a = {a}")
print(f"b = {b}")
print(f"c = {c}")
print("-" * 20)

# Step 3: Calculate the total energy E_total.
# E_total = integral from 0 to 1 of Φ(2x^2 + x) dx
# E_total = integral from 0 to 1 of [a*(2x^2+x)^2 + b*(2x^2+x) + c] dx
# After performing the integration symbolically, we get:
# E_total = a*(32/15) + b*(7/6) + c
term_a = a * Fraction(32, 15)
term_b = b * Fraction(7, 6)
term_c = c

E_total = term_a + term_b + term_c

# Step 4: Print the final equation and result.
print("The minimum total energy is calculated as:")
print(f"E_total = a*(32/15) + b*(7/6) + c")
print(f"E_total = ({a}) * (32/15) + ({b}) * (7/6) + ({c})")
print(f"E_total = ({term_a}) + ({term_b}) + ({term_c})")
print(f"E_total = {E_total}")

# The final answer in the required format
# print(f"\n<<<{E_total}>>>")