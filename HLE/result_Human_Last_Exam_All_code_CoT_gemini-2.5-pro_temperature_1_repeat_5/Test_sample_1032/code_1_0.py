import sys

# Rope length
L = 7 / 2

# 1. Area of the three unobstructed quadrants (Q1, Q2, Q4)
# Each is a triangle with area 1/2 * L^2
area_3_quads_num = 3 * (L**2)
area_3_quads_den = 2
area_3_quads = area_3_quads_num / area_3_quads_den

# 2. Area of the directly accessible region in Q3
# This is a 1x1 square
area_q3_direct = 1.0

# 3. Area accessible by wrapping around corners V1(-2,0) and V3(0,-2)
# Distance to these corners is 2.
dist_to_corner = 2
remaining_rope = L - dist_to_corner

# The area from each corner is a quarter of a taxi-cab diamond (a triangle).
# The area of a quarter-diamond is 1/2 * r^2
area_per_wrap = 0.5 * remaining_rope**2
area_q3_wrap = 2 * area_per_wrap # We have two such corners

# Convert to fractions for explanation
# L = 7/2
# area_3_quads = 3 * 1/2 * (7/2)^2 = 3 * 49/8 = 147/8
area_3_quads_num_f = 147
area_3_quads_den_f = 8

# remaining_rope = 7/2 - 2 = 3/2
# area_per_wrap = 1/2 * (3/2)^2 = 1/2 * 9/4 = 9/8
area_per_wrap_num_f = 9
area_per_wrap_den_f = 8

# area_q3_wrap = 2 * 9/8 = 18/8
area_q3_wrap_num_f = 18
area_q3_wrap_den_f = 8

# area_q3_direct = 1 = 8/8
area_q3_direct_num_f = 8
area_q3_direct_den_f = 8

# 4. Total Area
total_area_num_f = area_3_quads_num_f + area_q3_direct_num_f + area_q3_wrap_num_f
total_area_den_f = 8
total_area = total_area_num_f / total_area_den_f

# Print the final equation with all numbers
print("The total reachable area is the sum of several disjoint regions:")
print(f"1. Area in Q1, Q2, and Q4 = 3 * (1/2 * (7/2)^2) = {area_3_quads_num_f}/{area_3_quads_den_f}")
print(f"2. Directly reachable area in Q3 = 1 = {area_q3_direct_num_f}/{area_q3_direct_den_f}")
print(f"3. Area from wrapping around corner V1 = 1/2 * (7/2 - 2)^2 = {area_per_wrap_num_f}/{area_per_wrap_den_f}")
print(f"4. Area from wrapping around corner V3 = 1/2 * (7/2 - 2)^2 = {area_per_wrap_num_f}/{area_per_wrap_den_f}")
print("\nFinal Equation:")
print(f"Total Area = {area_3_quads_num_f}/{area_3_quads_den_f} + {area_q3_direct_num_f}/{area_q3_direct_den_f} + {area_per_wrap_num_f}/{area_per_wrap_den_f} + {area_per_wrap_num_f}/{area_per_wrap_den_f}")
print(f"Total Area = ({area_3_quads_num_f} + {area_q3_direct_num_f} + {area_q3_wrap_num_f}) / {total_area_den_f}")
print(f"Total Area = {total_area_num_f}/{total_area_den_f} = {total_area}")

sys.stdout.flush()
<<<21.625>>>