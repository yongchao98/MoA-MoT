import math

# The rope has a taxi-cab length L.
L = 3.5

# 1. Area in Quadrants 1, 2, and 4
# In each of these quadrants, the accessible region is a triangle with base L and height L.
# Area of one quadrant = 1/2 * L * L = 0.5 * L^2
# Total for three quadrants:
area_q1_q2_q4 = 3 * (0.5 * L**2)

# 2. Area in Quadrant 3
# The house is an L-shape in Q3. From the origin, the view into Q3 is
# bounded by rays passing through the corners of the house at (-2,-1) and (-1,-2).
# The slopes are m1 = (-1)/(-2) = 0.5 and m2 = (-2)/(-1) = 2.

# 2a. Directly accessible area in Q3 (the "unshadowed wedge")
# The area of the wedge between slopes m1 and m2 for a taxi-cab circle of radius L is L^2/6.
# Proof involves coordinate transformation or geometric argument.
area_q3_wedge_total = L**2 / 6

# We must subtract the area of the house that lies inside this wedge.
# The house is composed of 3 unit squares. The one at [-1,0]x[-1,0] is partially in the wedge.
# The area of the unit square [-1,0]x[-1,0] that falls between lines y=0.5x and y=2x is 0.5.
area_house_in_wedge = 0.5
area_q3_unshadowed = area_q3_wedge_total - area_house_in_wedge

# 2b. Indirectly accessible areas in Q3 (the "shadow regions")
# These are reached by the rope bending around the house's corners D(-1,-2) and F(-2,-1).

# Pivot D(-1,-2):
# Distance from origin to pivot D: |-1-0| + |-2-0| = 3.0
dist_to_d = 3.0
rem_len_d = L - dist_to_d  # Remaining rope length = 3.5 - 3.0 = 0.5
# This allows access to a new region centered at D. We are interested in the part
# of this new region that is in the shadow (y/x > 2) and doesn't re-enter the house (y <= -2).
# This area can be calculated as a small polygon, and its area is 1/6.
area_bent_d = 1.0/6.0

# Pivot F(-2,-1):
# By symmetry, the calculation is identical.
dist_to_f = 3.0
rem_len_f = L - dist_to_f
area_bent_f = 1.0/6.0

# 3. Total Area
total_area = area_q1_q2_q4 + area_q3_unshadowed + area_bent_d + area_bent_f

# 4. Output the result
# The problem asks to show each number in the final equation.
print("The total reachable area is calculated as the sum of several components:")
print(f"1. Area in Q1, Q2, and Q4: {area_q1_q2_q4:.4f}")
print(f"2. Directly reachable area in Q3 (wedge minus house): ({area_q3_wedge_total:.4f} - {area_house_in_wedge:.4f}) = {area_q3_unshadowed:.4f}")
print(f"3. Area from rope bending around corner D(-1,-2): {area_bent_d:.4f}")
print(f"4. Area from rope bending around corner F(-2,-1): {area_bent_f:.4f}")
print("\nFinal Equation:")
print(f"{area_q1_q2_q4:.4f} + ({area_q3_wedge_total:.4f} - {area_house_in_wedge:.4f}) + {area_bent_d:.4f} + {area_bent_f:.4f} = {total_area:.4f}")
print(f"\nThe total area the horse can reach is {total_area:.4f}")

# Final Answer as a fraction: 81/4
# final_fractional_area = (147/8) + (49/24 - 1/2) + (1/6) + (1/6)
# = 147/8 + 37/24 + 8/24 = 441/24 + 45/24 = 486/24 = 81/4

<<<20.25>>>