import math

# The user asks for the upper bound for the maximum norm (k_k,inf)
# in relation to the covolume (V).
# Based on interpretation, this relates to the systole of an arithmetic
# Riemann surface, which is bounded by its area (V).
# The inequality is: k_k,inf <= sqrt(2*V / sqrt(3)).
# This is equivalent to: k_k,inf <= (2^(1/2) / 3^(1/4)) * V^(1/2).

print("The upper bound is described by the inequality: k_k,inf <= C * V^p")
print("The final equation can be expressed as: k_k,inf <= (2^a / 3^b) * V^d")
print("\nThe numbers that define this equation are the bases (2, 3) and the exponents (a, b, d):")

a = 0.5
b = 0.25
d = 0.5

print(f"\n- Base 2 has exponent a = {a}")
print(f"- Base 3 has exponent b = {b} (in the denominator)")
print(f"- Base V (covolume) has exponent d = {d}")

# The problem asks to output each number in the final equation.
# The equation is k <= (2^0.5 / 3^0.25) * V^0.5
# So the numbers are 2, 0.5, 3, 0.25, 0.5
print("\nExplicitly, the numbers in the final equation k_k,inf <= (2^0.5 / 3^0.25) * V^0.5 are:")
print("Base: 2")
print("Exponent: 0.5")
print("Base: 3")
print("Exponent: 0.25")
print("Base: V")
print("Exponent: 0.5")
