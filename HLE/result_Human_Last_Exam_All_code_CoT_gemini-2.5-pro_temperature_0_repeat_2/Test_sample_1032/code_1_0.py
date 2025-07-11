import math

# 1. Define constants
rope_length = 7.0 / 2.0

# 2. Calculate the total area of the reachable diamond shape without obstacles
# The region is |x| + |y| <= L, which is a square with vertices (L,0), (0,L), (-L,0), (0,-L).
# Its area is 2 * L^2.
total_diamond_area = 2 * rope_length**2

# 3. Define the area of the house obstacle
# The house consists of three unit squares.
house_area = 3.0

# 4. Calculate the area lost due to the "shadow" effect of the house corners.
# The rope must bend around the outer corners C1=(-2,-1) and C2=(-1,-2).
# Let's calculate the area reduction for the corner C2=(-1,-2).
# This affects the region where x is in [-1, 0] and y <= -2.

# First, calculate the area in this shadow zone if the path were standard (straight).
# This is the area of {(x,y) | -1<=x<=0, y<=-2, and -x-y<=L}.
# This forms a trapezoid with vertices (0, -3.5), (0, -2), (-1, -2), (-1, -2.5).
# Its area is 0.5 * (side1 + side2) * width = 0.5 * ((3.5-2) + (2.5-2)) * 1 = 1.0.
area_lost_standard_path = 1.0

# Second, calculate the actual reachable area with the bent rope.
# The path length is d(O,C2) + d(C2,P) = (|-1|+|-2|) + (|x+1|+|y+2|) = 3 + (x+1) - (y+2) = x-y+2.
# The condition is x-y+2 <= 3.5, or y >= x-1.5.
# This forms a triangle with vertices (-1, -2.5), (-1, -2), (-0.5, -2).
# Its area is 0.5 * base * height = 0.5 * 0.5 * 0.5 = 0.125.
area_lost_actual_path = 0.125

# The reduction in area due to this one corner is the difference.
area_reduction_per_corner = area_lost_standard_path - area_lost_actual_path

# The total reduction is for two corners (C1 and C2).
total_area_reduction = 2 * area_reduction_per_corner

# 5. Calculate the final area
# Final Area = Total Diamond Area - House Area - Total Shadow Area Reduction
final_area = total_diamond_area - house_area - total_area_reduction

# 6. Print the final calculation as a single equation
print("The total reachable area is calculated by taking the area of the full taxi-cab circle, subtracting the area of the house, and subtracting the area lost in the 'shadows' cast by the house's corners.")
print("Final Equation:")
print(f"{total_diamond_area} - {house_area} - 2 * ({area_lost_standard_path} - {area_lost_actual_path}) = {final_area}")

<<<19.75>>>