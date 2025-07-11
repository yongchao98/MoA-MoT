# Rope's taxi-cab length
L = 7 / 2

# 1. Area in Quadrants 1, 2, and 4 (unobstructed)
# The area in each quadrant is a triangle with area 0.5 * L^2
area_per_unobstructed_quadrant = 0.5 * L**2
area_q1_q2_q4 = 3 * area_per_unobstructed_quadrant

# 2. Area in Quadrant 3 (obstructed)
# We split Q3 by the line y=x.

# 2a. Sub-region where y > x
# The total possible area is a triangle with vertices (0,0), (0,-L), (-L/2,-L/2).
# Its area is L^2 / 4.
total_area_y_gt_x = L**2 / 4
# The portion of the house in this region is the sum of:
# - Half of the square [-1,0]x[-1,0] (Area = 0.5)
# - All of the square [-2,-1]x[-1,0] (Area = 1.0)
house_area_y_gt_x = 0.5 + 1.0
reachable_area_y_gt_x = total_area_y_gt_x - house_area_y_gt_x

# 2b. Sub-region where x > y (symmetric)
total_area_x_gt_y = L**2 / 4
# The portion of the house in this region is the sum of:
# - Half of the square [-1,0]x[-1,0] (Area = 0.5)
# - All of the square [-1,0]x[-2,-1] (Area = 1.0)
house_area_x_gt_y = 0.5 + 1.0
reachable_area_x_gt_y = total_area_x_gt_y - house_area_x_gt_y

# Total area in Q3 is the sum of the two sub-regions
area_q3 = reachable_area_y_gt_x + reachable_area_x_gt_y

# 3. Total Reachable Area
total_area = area_q1_q2_q4 + area_q3

# --- Output the final equation and result ---
print("The total reachable area is the sum of the areas from all four quadrants.")
print(f"Area(Q1+Q2+Q4) = 3 * (0.5 * {L}**2) = {area_q1_q2_q4}")
print(f"Area(Q3) = ({total_area_y_gt_x} - {house_area_y_gt_x}) + ({total_area_x_gt_y} - {house_area_x_gt_y}) = {area_q3}")
print(f"Total Area = {area_q1_q2_q4} + {area_q3} = {total_area}")

<<<21.5>>>