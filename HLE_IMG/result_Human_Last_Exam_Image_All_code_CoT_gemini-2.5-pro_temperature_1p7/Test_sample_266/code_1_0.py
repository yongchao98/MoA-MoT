import math

# Define the side length of the red hexagon
s = 3

# The formula for the area of a regular hexagon is (3 * sqrt(3) / 2) * s^2
# We will calculate this value.
area_hexagon = (3 * math.sqrt(3) / 2) * s**2

print("The edge length of the red hexagon is given as s = 3.")
print("The area of a regular hexagon is calculated using the formula: (3 * sqrt(3) / 2) * s^2")
print(f"Plugging in the value of s = {s}:")
print(f"Area = (3 * sqrt(3) / 2) * {s}^2")
print(f"Area = (3 * {math.sqrt(3):.4f} / 2) * {s**2}")
print(f"Area = {area_hexagon:.2f}")