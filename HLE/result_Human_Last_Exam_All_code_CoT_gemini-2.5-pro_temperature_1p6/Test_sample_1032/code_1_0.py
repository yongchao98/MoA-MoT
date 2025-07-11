from fractions import Fraction

# The total taxi-cab length of the rope
R = Fraction(7, 2)

# 1. Calculate the area in Quadrants 1, 2, and 4.
# Each is a right triangle with legs of length R.
# The area of one such triangle is 1/2 * R * R.
area_one_quadrant = Fraction(1, 2) * R**2
area_q1_q2_q4 = 3 * area_one_quadrant

# 2. Calculate the reachable area in Quadrant 3.
# This area is composed of several parts.

# 2a. Directly reachable area next to the y-axis.
# This is a triangle bounded by the y-axis, the main boundary |x|+|y|=R,
# and the line from the origin to corner C2=(-1,-2), which is y=2x.
# The base on the y-axis is R. The height is the absolute x-coordinate
# where y=2x intersects -x-y=R, which is |-R/3| = R/3.
area_q3_a = Fraction(1, 2) * R * (R / 3)

# 2b. Directly reachable area next to the x-axis.
# This is a symmetric triangle bounded by the x-axis, the main boundary,
# and the line from the origin to corner C1=(-2,-1), which is y=x/2.
# The base on the x-axis is R. The height is the absolute y-coordinate
# where y=x/2 intersects -x-y=R, which is |-R/3| = R/3.
area_q3_b = Fraction(1, 2) * R * (R / 3)

# 2c. Area reached by wrapping the rope around corner C1=(-2,-1).
# The taxi-cab distance to C1 is |-2|+|-1|=3.
dist_to_c1 = Fraction(3)
remaining_rope_c1 = R - dist_to_c1
# This forms a new pivot. The reachable area is a quarter of a diamond
# with radius equal to the remaining rope length. Total diamond area is 2*r^2.
area_wrap_c1 = Fraction(1, 4) * (2 * remaining_rope_c1**2)

# 2d. Area reached by wrapping the rope around corner C2=(-1,-2).
# The taxi-cab distance to C2 is |-1|+|-2|=3.
dist_to_c2 = Fraction(3)
remaining_rope_c2 = R - dist_to_c2
# This area is symmetric to the C1 case.
area_wrap_c2 = Fraction(1, 4) * (2 * remaining_rope_c2**2)

# 3. Sum all the parts for the total area.
total_area = area_q1_q2_q4 + area_q3_a + area_q3_b + area_wrap_c1 + area_wrap_c2

# 4. Print the final equation.
print(f"Total Area = (Area Q1+Q2+Q4) + (Area Q3 direct a) + (Area Q3 direct b) + (Area Q3 wrapped C1) + (Area Q3 wrapped C2)")
print(f"Total Area = {area_q1_q2_q4} + {area_q3_a} + {area_q3_b} + {area_wrap_c1} + {area_wrap_c2}")
print(f"Total Area = {float(area_q1_q2_q4):.4f} + {float(area_q3_a):.4f} + {float(area_q3_b):.4f} + {float(area_wrap_c1):.4f} + {float(area_wrap_c2):.4f}")
print(f"Final Equation: {area_q1_q2_q4} + {area_q3_a} + {area_q3_b} + {area_wrap_c1} + {area_wrap_c2} = {total_area}")
print(f"The total area is {total_area} square units, which is approximately {float(total_area):.4f}.")
