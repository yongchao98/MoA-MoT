import numpy as np

# Coordinates of the five leg tips in a plane
P1 = (0, 0)
P2 = (2, 0)
P3 = (2, 2)
P4 = (0, 2)
P5 = (1, 4)

# The first four points P1, P2, P3, P4 form a square.
# We find the equation of the circle that passes through these four points.
# The center of the circle is the geometric center of the square.
center_x = (P1[0] + P3[0]) / 2.0
center_y = (P1[1] + P3[1]) / 2.0

# The radius squared is the squared distance from the center to any of the four corner points.
# Using P1=(0,0) and center=(1.0, 1.0)
radius_sq = (P1[0] - center_x)**2 + (P1[1] - center_y)**2

print("Step 1: Determine the circle for the four rectangular legs.")
print(f"The circle passing through P1, P2, P3, and P4 has:")
print(f"Equation: (x - {center_x})^2 + (y - {center_y})^2 = {radius_sq}")
print("-" * 50)

print("Step 2: Check if the fifth leg P5(1,4) lies on this circle.")
x5, y5 = P5
# Calculate the left-hand side (LHS) of the circle equation for P5
lhs_p5 = (x5 - center_x)**2 + (y5 - center_y)**2

print(f"Substitute P5's coordinates into the equation's left side:")
print(f"({x5} - {center_x})^2 + ({y5} - {center_y})^2 = {lhs_p5}")
print(f"We compare this result ({lhs_p5}) to the radius squared ({radius_sq}).")

if np.isclose(lhs_p5, radius_sq):
    print("The values are equal. The five points are concyclic.")
else:
    print("The values are not equal. The five points are NOT concyclic.")

print("-" * 50)
print("Conclusion:")
print("The geometric analysis shows the five leg positions are not concyclic.")
print("For a chair with non-concyclic legs, mathematical theorems state:")
print("1. It is IMPOSSIBLE to place it on a PERFECT sphere (0 solutions).")
print("2. It is ALWAYS POSSIBLE to place it on a 'smooth but uneven' sphere in at least TWO locations.")
print("Since the problem specifies an uneven surface, the minimum cardinality of the set of locations is 2.")
