from fractions import Fraction

# The problem asks for the area of a region defined by a taxi-cab distance
# with an obstacle. We calculate this by summing the areas of disjoint regions
# that become accessible as the rope moves around the obstacle.

# Step 1: Define the rope length
L = Fraction(7, 2)

# Step 2: Calculate the area in the three unobstructed quadrants (Q1, Q2, Q4).
# This is 3/4 of the total diamond area (2 * L^2).
area_main = 3 * (Fraction(1, 2) * L**2)

# Step 3: Calculate the area in Q3 reachable by a straight, unblocked path.
# The view from the origin is blocked by house corners, leaving two accessible
# triangular sectors. The total area of these two sectors is L^2 / 3.
area_q3_straight = (L**2) / 3

# Step 4: Calculate the area accessible by the rope bending around the outer corners.
# The path to corners (-2,0) or (0,-2) has a taxi-cab length of 2.
d_outer = Fraction(2)
L_rem = L - d_outer # Remaining rope length is 1.5

# For corner (-2,0), a half-diamond of area (L_rem)^2 is added.
area_c1 = L_rem**2
# For corner (0,-2), another half-diamond of area (L_rem)^2 is added.
area_c2 = L_rem**2

# Step 5: Calculate the area accessible by the rope bending into the concave corner (-1,-1).
# The shortest path along the house boundary to (-1,-1) is also 2.
# This opens up a quarter-diamond (a triangle) in the nook of the house.
area_c3 = Fraction(1, 2) * L_rem**2

# Step 6: Sum the disjoint areas for the total area.
total_area = area_main + area_q3_straight + area_c1 + area_c2 + area_c3

# Step 7: Print the final calculation as requested.
# The final equation shows each calculated component area.
print("The total reachable area is the sum of five disjoint regions:")
print(f"Area = (Area in Q1,2,4) + (Area in Q3, straight) + (Area around (-2,0)) + (Area around (0,-2)) + (Area in nook (-1,-1))")
print(f"Total Area = {area_main.numerator}/{area_main.denominator} + {area_q3_straight.numerator}/{area_q3_straight.denominator} + {area_c1.numerator}/{area_c1.denominator} + {area_c2.numerator}/{area_c2.denominator} + {area_c3.numerator}/{area_c3.denominator} = {total_area.numerator}/{total_area.denominator}")
print(f"\nThe final area is {total_area.numerator}/{total_area.denominator}, which is approximately {float(total_area):.4f} square units.")
