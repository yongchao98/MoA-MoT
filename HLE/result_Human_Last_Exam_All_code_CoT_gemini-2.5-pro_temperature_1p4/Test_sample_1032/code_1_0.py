import math

# The rope has a taxi-cab length L.
L = 7 / 2

# --- Step 1: Calculate the area in the three unobstructed quadrants (Q1, Q2, Q4) ---
# The area of the full diamond |x|+|y|<=L is 2*L^2.
# The area in one quadrant is 0.5 * L^2.
# We have three such quadrants.
area_unobstructed = 3 * (0.5 * L**2)

# --- Step 2: Calculate the area in Quadrant 3, which is obstructed ---
# The area is split into three parts based on the pivots the rope must bend around.

# --- Part 2a: Area reachable by bending around pivot P1=(-2, 0) ---
# This applies to the region where x <= -2.
# The taxi-cab distance from the origin (0,0) to P1 is |-2-0| + |0-0| = 2.
dist_to_p1 = 2.0
rope_left_p1 = L - dist_to_p1
# The reachable area from P1 is a diamond |x+2| + |y| <= rope_left_p1.
# We are interested in the quarter of this diamond in the third quadrant (x<=-2, y<=0).
# The area of this quarter-diamond is 0.5 * rope_left_p1^2.
area_p1_region = 0.5 * rope_left_p1**2

# --- Part 2b: Area reachable by bending around pivot P2=(0, -2) ---
# This is symmetric to the P1 case and applies to y <= -2.
# The distance to P2 is |0-0| + |-2-0| = 2.
dist_to_p2 = 2.0
rope_left_p2 = L - dist_to_p2
# The area of the corresponding quarter-diamond is 0.5 * rope_left_p2^2.
area_p2_region = 0.5 * rope_left_p2**2

# --- Part 2c: Area reachable in the square [-2, -1] x [-2, -1] ---
# The path to this region bends around the concave corner C=(-1, -1).
# The distance to C is |-1-0| + |-1-0| = 2.
dist_to_c = 2.0
rope_left_c = L - dist_to_c
# The reachable points must satisfy |x+1| + |y+1| <= rope_left_c inside the square.
# The square has an area of 1.
# A small triangle at the corner (-2, -2) is unreachable.
# The line |x+1|+|y+1| = 1.5, or x+y+2 = -1.5 => x+y = -3.5, defines the boundary.
# Re-evaluating in terms of coordinates relative to C (-1, -1):
# x' = x+1, y' = y+1. The region is [-1, 0]x[-1, 0].
# The condition is |x'|+|y'| <= 1.5. The excluded region is |x'|+|y'| > 1.5.
# In the square's domain, this means -x' -y' > 1.5 or x'+y' < -1.5.
# This cuts off a triangle with vertices (-1,-1), (-1,-0.5), (-0.5,-1) in (x',y') coords.
# The legs of this triangle have length 0.5.
area_cutoff_triangle = 0.5 * (1.5 - 1)**2 # A more general way to see it is 0.5*(rope_left_c - corner_dist_from_c)^2. Here dist is 1.
area_cutoff_triangle = 0.5 * 0.5 * 0.5
area_c_region = 1.0 - area_cutoff_triangle

# --- Step 3: Sum the areas and print the results ---
total_area = area_unobstructed + area_p1_region + area_p2_region + area_c_region

print("The total reachable area is the sum of four distinct regions:")
print(f"1. Area in the three unobstructed quadrants (Q1, Q2, Q4): {area_unobstructed}")
print(f"2. Area in Q3 reachable via pivot P1=(-2,0): {area_p1_region}")
print(f"3. Area in Q3 reachable via pivot P2=(0,-2): {area_p2_region}")
print(f"4. Area in Q3 reachable via pivot C=(-1,-1): {area_c_region}")
print("-" * 20)
print(f"Total Area = {area_unobstructed} + {area_p1_region} + {area_p2_region} + {area_c_region}")
print(f"Total Area = {total_area}")

<<<21.5>>>