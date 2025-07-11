import math

# The problem reduces to solving for r_0 in the equation:
# (3 * r_0 - 37) / (r_0 + 4) = 1 / sqrt(2)
# The exact solution for r_0 is of the form (a + b * sqrt(c)) / d.
a = 226
b = 49
c = 2
d = 17

# Calculate the numerical value of r_0
sqrt_c = math.sqrt(c)
r_0_numerical = (a + b * sqrt_c) / d

# Print the final equation with its components
print("The radial distance r_0 is determined by the equation:")
print(f"r_0 = ({a} + {b} * sqrt({c})) / {d}")
print("\nThe numbers in the final equation are:")
print(f"a = {a}")
print(f"b = {b}")
print(f"c = {c}")
print(f"d = {d}")

# Print the final numerical answer for r_0
print(f"\nThe numerical value of the radial distance r_0 is:")
print(r_0_numerical)