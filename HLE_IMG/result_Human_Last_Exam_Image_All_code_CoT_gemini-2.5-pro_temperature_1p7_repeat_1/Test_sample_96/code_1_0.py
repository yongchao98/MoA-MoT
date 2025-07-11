import math

# Define constants
pi = math.pi
sqrt3 = math.sqrt(3)

# 1. Area of the circle with radius 2
r = 2
area_circle = pi * r**2
print(f"The area of the circle x^2 + y^2 = 4 is π * {r}^2 = {4}π.")
print("-" * 30)

# 2. Calculate the area of the region R inside the circle
# This region is shaded gray in the provided image and consists of three parts.

# Part 1: Intersection with the rectangle R_2 = [-1, 0] x [0, 1]
# This rectangle is fully inside the circle.
area_part1 = 1.0
print("The area of the first part of the intersection (the rectangle [-1, 0] x [0, 1]) is 1.")

# Part 2: Intersection with the rectangle R_1 = [-2, -1] x [0, 2]
# This is the area under the circle y=sqrt(4-x^2) from x=-2 to x=-1.
# The integral evaluates to 2π/3 - sqrt(3)/2.
area_part2 = (2 * pi / 3) - (sqrt3 / 2)
print(f"The area of the second part of the intersection is (2π/3 - √3/2).")

# Part 3: Intersection with the rectangle R_3 = [1, 2] x [-1, 0]
# This results in a smaller rectangle [1, √3] x [-1, 0].
area_part3 = sqrt3 - 1
print(f"The area of the third part of the intersection is (√3 - 1).")
print("-" * 30)

# 3. Sum the intersection areas
total_intersection_area = area_part1 + area_part2 + area_part3
print("The total area of the intersection (region R inside the circle) is the sum of these parts:")
print(f"Area_intersection = {area_part1:.4f} + ({2/3:.4f}π - {sqrt3/2:.4f}) + ({sqrt3:.4f} - 1)")
print(f"Area_intersection = 1 + {area_part2:.4f} + {area_part3:.4f}")
print(f"Area_intersection = {total_intersection_area:.4f}")
print("Symbolically, Area_intersection = 1 + (2π/3 - √3/2) + (√3 - 1) = 2π/3 + √3/2.")
print("-" * 30)

# 4. Calculate the final area (Area of circle - Area of intersection)
final_area = area_circle - total_intersection_area
print("The area inside the circle but outside of region R is:")
print("Area = Area_circle - Area_intersection")
print(f"Area = 4π - (2π/3 + √3/2)")
print(f"Area = 12π/3 - 2π/3 - √3/2")
print(f"Area = 10π/3 - √3/2")
print("-" * 30)
print(f"Final calculated area: {final_area:.6f}")
print("The final equation is:")
print(f"{4} * π - ({1} + ({2} * π / {3} - {3}**0.5 / {2}) + ({3}**0.5 - {1})) = {10} * π / {3} - {3}**0.5 / {2}")
<<<10*pi/3 - sqrt(3)/2>>>