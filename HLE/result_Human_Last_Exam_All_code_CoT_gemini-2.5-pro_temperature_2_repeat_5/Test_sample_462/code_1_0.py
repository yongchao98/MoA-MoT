import math

# Define the side length of the square
side_length = 20

# Step 1: Determine the radius of the circles.
# The diameters of two circles span the side of the square.
# 2 * diameter = side_length => 2 * (2 * radius) = side_length => 4 * radius = side_length
radius = side_length / 4.0

# Step 2: Calculate the area of the square.
area_square = side_length ** 2

# Step 3: Calculate the total area of the four circles.
area_four_circles = 4 * math.pi * (radius ** 2)

# Step 4: Calculate the area of the region between the circles.
area_between = area_square - area_four_circles

# Step 5: Print the steps and the final rounded answer.
print(f"The side length of the square is {side_length} cm.")
print(f"The radius of each of the four circles is {side_length} / 4 = {radius} cm.")
print(f"First, we calculate the area of the square: {side_length} * {side_length} = {area_square} cm^2.")
print(f"Next, we calculate the total area of the four circles: 4 * pi * {radius}^2 = {area_four_circles} cm^2.")
print("Finally, we subtract the area of the circles from the area of the square.")
print(f"The final equation is: Area = {area_square} - {area_four_circles}")
print(f"The area of the region between the circles is {area_between} cm^2.")
print(f"Rounded to the nearest hundredth, the area is {round(area_between, 2)} cm^2.")
