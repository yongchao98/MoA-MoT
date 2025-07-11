import math

# From the derivation, the values (x_0, y_0, z_0) form a sphere
# with radius squared R^2 = 0.5 * 10^25.
r_squared_mantissa = 0.5
r_squared_exponent = 25
R_squared = r_squared_mantissa * (10**r_squared_exponent)

# The surface area of a sphere is given by the formula A = 4 * pi * R^2.
four_pi = 4 * math.pi
area = four_pi * R_squared

# Output the calculation steps as requested.
print("The values (x_0, y_0, z_0) satisfy the equation of a sphere: x_0^2 + y_0^2 + z_0^2 = R^2")
print(f"where R^2 = {r_squared_mantissa} * 10^{r_squared_exponent} = {R_squared}")
print("\nThe area of this sphere is calculated using the formula: A = 4 * pi * R^2")
print(f"A = 4 * {math.pi} * {R_squared}")
print(f"A = {four_pi} * {R_squared}")
print(f"\nThe final area is: {area}")
