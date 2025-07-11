import math

# 1. Define constants
pi = math.pi
sqrt3 = math.sqrt(3)
r = 2

# 2. Calculate the total area of the circle
area_circle = pi * r**2
print(f"The total area of the circle with radius {r} is: {area_circle:.4f} (which is 4π)")

# 3. Calculate the area of the intersection of region R and the circle
# This area is the sum of three parts as explained in the derivation.
# Part 1: A 1x1 rectangle in the second quadrant.
area_part1 = 1.0
# Part 2: The area under the circle's arc from x=-2 to x=-1.
area_part2 = (2 * pi / 3) - (sqrt3 / 2)
# Part 3: The area of R inside the circle in the fourth quadrant.
area_part3 = (pi / 3) + (sqrt3 / 2) - 1.0

# Total intersection area
area_intersection = area_part1 + area_part2 + area_part3
print(f"The area of intersection (R inside the circle) is the sum of three parts:")
print(f"  Part 1 (Rectangle): {area_part1:.4f}")
print(f"  Part 2 (Q2 Arc): {area_part2:.4f}")
print(f"  Part 3 (Q4 Region): {area_part3:.4f}")
print(f"Total Intersection Area = {area_intersection:.4f} (which is exactly π)")

# 4. Calculate the final area (inside the circle, outside of R)
final_area = area_circle - area_intersection

# 5. Print the final equation and result
# The problem asks to output each number in the final equation.
print("\nThe final calculation is:")
print(f"Area = (Area of Circle) - (Area of Intersection)")
print(f"Area = 4 * π - π")
print(f"Area = 3 * π")
print(f"The numerical value is approximately: {final_area:.4f}")

# Final answer in the required format
# The exact answer is 3π.
# final_answer_value = 3 * pi
# print(f'<<<{final_answer_value}>>>')