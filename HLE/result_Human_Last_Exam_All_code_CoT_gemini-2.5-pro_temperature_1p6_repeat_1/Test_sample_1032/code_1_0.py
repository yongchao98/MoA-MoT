# The taxi-cab length of the rope
L = 7.0 / 2.0

# 1. Calculate the area in Quadrants 1, 2, and 4.
# The area of the full diamond is 2 * L^2.
# The area in one quadrant is a triangle with area L^2 / 2.
area_q1 = (L**2) / 2.0
area_q2 = (L**2) / 2.0
area_q4 = (L**2) / 2.0
area_unobstructed_quadrants = area_q1 + area_q2 + area_q4

# 2. Calculate the area in Quadrant 3.
# This is composed of a directly reachable part and parts reachable by wrapping around corners.

# 2a. Directly reachable area in Q3.
# This is the area in { (x,y) | |x|+|y|<=L, x<=0, y<=0, (x>=-1 or y>=-1) } minus the house area.
# Area for x >= -1 strip: The region is a trapezoid with vertices (-1,0), (0,0), (0,-3.5), (-1,-2.5).
# Area of this trapezoid = 0.5 * (2.5 + 3.5) * 1 = 3.
# The house occupies two 1x1 squares in this strip. Area = 2.
# Net direct area in this strip = 3 - 2 = 1.
area_direct_strip_x = 1.0

# Area for y >= -1 strip is symmetric. Net direct area = 1.
area_direct_strip_y = 1.0

# The overlap is the square [-1,0]x[-1,0]. Area = 1. The house occupies this entire square. So net overlap is 0.
area_direct_overlap = 0.0

# Total directly reachable area in Q3 using inclusion-exclusion.
area_direct_q3 = area_direct_strip_x + area_direct_strip_y - area_direct_overlap

# 2b. Area reachable by wrapping around convex corners.
# For a corner V, path length to it is d, remaining rope is r = L-d.
# The new area is a quarter-diamond with area r^2 / 2.

# Corner A(-2,0):
d_A = 2.0
r_A = L - d_A
area_A = (r_A**2) / 2.0

# Corner C(0,-2):
d_C = 2.0
r_C = L - d_C
area_C = (r_C**2) / 2.0

# Corner F(-2,-1): Path is O->A->F, d = d(O,A)+d(A,F) = 2+1=3.
d_F = 3.0
r_F = L - d_F
area_F = (r_F**2) / 2.0

# Corner D(-1,-2): Path is O->C->D, d = d(O,C)+d(C,D) = 2+1=3.
d_D = 3.0
r_D = L - d_D
area_D = (r_D**2) / 2.0

# 3. Sum all the areas.
total_area = area_unobstructed_quadrants + area_direct_q3 + area_A + area_C + area_F + area_D

# Print the calculation steps
print("Calculation of the reachable area:")
print(f"Rope length (L): {L}")
print("-" * 30)
print("1. Area in unobstucted Quadrants (Q1, Q2, Q4):")
print(f"   Area = 3 * (L^2 / 2) = 3 * ({L}^2 / 2) = {area_unobstructed_quadrants}")
print("-" * 30)
print("2. Area in obstructed Quadrant (Q3):")
print(f"   a. Directly reachable area = {area_direct_q3}")
print(f"   b. Area around corner A(-2,0) = (L - 2)^2 / 2 = ({L-d_A:.1f})^2 / 2 = {area_A}")
print(f"   c. Area around corner C(0,-2) = (L - 2)^2 / 2 = ({L-d_C:.1f})^2 / 2 = {area_C}")
print(f"   d. Area around corner F(-2,-1) = (L - 3)^2 / 2 = ({L-d_F:.1f})^2 / 2 = {area_F}")
print(f"   e. Area around corner D(-1,-2) = (L - 3)^2 / 2 = ({L-d_D:.1f})^2 / 2 = {area_D}")
print("-" * 30)
print("3. Total Area:")
print("   Total Area = Area(Q1,Q2,Q4) + Area(Direct Q3) + Area(A) + Area(C) + Area(F) + Area(D)")
print(f"   Total Area = {area_unobstructed_quadrants} + {area_direct_q3} + {area_A} + {area_C} + {area_F} + {area_D}")
print(f"   Total Area = {total_area}")

<<<22.875>>>