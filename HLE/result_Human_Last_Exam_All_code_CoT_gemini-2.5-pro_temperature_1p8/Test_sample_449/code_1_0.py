import math

# Starting point
x1 = 3000
x2 = 4000

# Magnitude of the starting point vector
r = math.sqrt(x1**2 + x2**2)

# cos(theta) and sin(theta)
cos_theta = x1 / r
sin_theta = x2 / r

# cos(2*theta)
cos_2_theta = cos_theta**2 - sin_theta**2

# Probability approximation using the formula: -1/(pi*r^4) * cos(2*theta)
# This formula is for P(T_0 < T_A), which is the quantity we need.
prob = -1 / (math.pi * r**4) * cos_2_theta

# We want the result with two significant digits.
print(f"The starting point is ({x1},{x2}).")
print(f"The distance from the origin is r = {r:.0f}.")
print(f"cos(theta) = {cos_theta}, sin(theta) = {sin_theta}.")
print(f"cos(2*theta) = {cos_2_theta:.4f}.")
print(f"The probability is approximately {prob:.2g}.")
print(f"\nFinal calculation in full:")
print(f"P = - (1 / (pi * {r:.0f}^4)) * ({cos_2_theta:.4f})")
print(f"P = (1 / (pi * {r**4:.4g})) * ({-cos_2_theta:.4f})")
print(f"P = {prob:.4g}")