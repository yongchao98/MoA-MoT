import math

# Rope length
L = 7.0 / 2.0

# 1. Area in Quadrants 1, 2, and 4
# Each quadrant is a triangle with base L and height L. Area = 0.5 * L * L
area_q1_q2_q4 = 3 * (0.5 * L * L)

# 2. Area in Quadrant 3, partitioned into regions
# The house itself is unreachable.

# 2a. Region S4: The square [-2, -1] x [-2, -1]
# The condition for being reachable is |x-y| >= 0.5.
# We calculate the area of the excluded region |x-y| < 0.5 within the unit square.
# This excluded region consists of two identical triangles.
# Each triangle has a base of 0.5 and height of 0.5.
area_excluded_in_s4 = 2 * (0.5 * 0.5 * 0.5)
area_s4_reachable = 1 - area_excluded_in_s4

# 2b. Region R_N: x <= -2, y in [-1, 0]
# Reachable area is where -x-y <= 3.5. This forms a trapezoid.
# Vertices: (-2,0), (-3.5,0), (-2.5,-1), (-2,-1).
# Area = 0.5 * (base1 + base2) * height = 0.5 * (1.5 + 0.5) * 1 = 1.0
area_rn_reachable = 1.0

# 2c. Region R_W: y <= -2, x in [-1, 0]
# By symmetry with R_N, the area is the same.
area_rw_reachable = 1.0

# 2d. Region R_C: x <= -2, y <= -2
# Reachable area is where -x-y <= 3.5. This forms a triangle.
# Vertices: (-2, -1.5), (-1.5, -2), (-2, -2)
# Area = 0.5 * base * height = 0.5 * 0.5 * 0.5
area_rc_reachable = 0.5 * 0.5 * 0.5

# Total area in Q3
area_q3 = area_s4_reachable + area_rn_reachable + area_rw_reachable + area_rc_reachable

# Total reachable area
total_area = area_q1_q2_q4 + area_q3

# Output the equation
print("The total area is calculated as the sum of areas in each region:")
print("Area = (Area in Q1, Q2, Q4) + (Reachable Area in Q3)")
print("Area in Q1, Q2, Q4 = 3 * (1/2 * (7/2)^2) = {}".format(area_q1_q2_q4))
print("Reachable Area in Q3 is the sum of areas in its sub-regions:")
print("  - Area in square [-2,-1]x[-2,-1] = 1 - 2*(1/2 * 1/2 * 1/2) = {}".format(area_s4_reachable))
print("  - Area in region x<=-2, y in [-1,0] = {}".format(area_rn_reachable))
print("  - Area in region y<=-2, x in [-1,0] = {}".format(area_rw_reachable))
print("  - Area in region x<=-2, y<=-2 = 1/2 * 1/2 * 1/2 = {}".format(area_rc_reachable))
print("Total Reachable Area in Q3 = {} + {} + {} + {} = {}".format(
    area_s4_reachable, area_rn_reachable, area_rw_reachable, area_rc_reachable, area_q3))
print("\nFinal Equation:")
print("Total Area = 3 * (1/2 * (7/2)^2) + (1 - 2 * (1/2 * 1/2 * 1/2)) + 1.0 + 1.0 + (1/2 * 1/2 * 1/2)")
print("Total Area = {} + {} + {} + {} + {} = {}".format(
    area_q1_q2_q4, area_s4_reachable, area_rn_reachable, area_rw_reachable, area_rc_reachable, total_area))
print("\nFinal Answer:")
print(total_area)