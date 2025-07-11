import math

# This problem requires finding the optimal shape of a 1 m^3 volume of uniform
# density to maximize the gravitational field at a given point A.
# The solution to this problem is a specific shape whose boundary is defined
# in spherical coordinates (with A at the origin) by the equation r^2 = cos(θ)/C,
# where C is a constant determined by the volume constraint.

# The derivation involves the following steps:
# 1. The volume V of this shape is calculated as: V = 4 * math.pi / (15 * C**(3/2)).
# 2. We are given V = 1 m^3, so we can determine C.
# 3. The furthest point on the surface from A, r_max, occurs at θ=0.
#    The distance is given by r_max = 1 / sqrt(C).
# 4. By combining these equations, we derive the final formula for the furthest distance:
#    r_max = (15 / (4 * π))^(1/3).

# Here, we define the numbers present in this final equation.
numerator = 15
denominator_constant = 4
power = 1/3

# We use the high-precision value of π from the math library.
pi_value = math.pi

print(f"The calculation for the furthest distance (r_max) is based on the formula:")
print(f"r_max = ({numerator} / ({denominator_constant} * π)) ** ({power})")
print("-" * 30)

# Perform the calculation
result = (numerator / (denominator_constant * pi_value)) ** power

# Print the final result
print("The calculated furthest distance is:")
print(result)