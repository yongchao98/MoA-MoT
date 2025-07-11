import math

# This script calculates the required ratio in its exact, symbolic form.
# Let s be the side length of the cube.

print("This program calculates the ratio of the specified areas on a cube's surface.")
print("The final answer will be presented in its exact form, including terms like pi and sqrt(3).")
print("-" * 30)

# Step 1: Total Surface Area (Area_S)
# A cube has 6 faces, each of area s^2.
# Area_S = 6 * s^2
print("Step 1: Total Surface Area")
print("The total surface area of a cube with side length s is Area(S) = 6 * s^2.")
print("-" * 30)

# Step 2: Area on the 3 faces adjacent to vertex P (Area_adj)
# For any point on a face adjacent to the vertex P, the shortest surface distance
# is the straight-line distance on that face. The maximum such distance is to the
# corner opposite P on that face, which is sqrt(s^2 + s^2) = sqrt(2)*s.
# Since the distance condition is d <= sqrt(2)*s, the entire area of these three
# adjacent faces is included in the region D.
# Area_adj = 3 * s^2
print("Step 2: Area on Adjacent Faces")
print("The entire area of the 3 faces adjacent to vertex P is included in region D.")
print("The contribution from these faces is Area_adj = 3 * s^2.")
print("-" * 30)

# Step 3: Area on one of the 3 faces opposite to vertex P (A_opp)
# For a point on an opposite face, the shortest path requires unfolding the cube.
# The area of D on one such face can be found through integration, which yields:
# A_opp = s^2 * (pi/3 - 3/2 + sqrt(3)/2)
# The total contribution from the three opposite faces is 3 * A_opp.
print("Step 3: Area on Opposite Faces")
print("The area on one opposite face, A_opp, is calculated to be s^2 * (pi/3 - 3/2 + sqrt(3)/2).")
print("The total contribution from the 3 opposite faces is 3 * A_opp = s^2 * (pi - 9/2 + 3*sqrt(3)/2).")
print("-" * 30)

# Step 4: Total Area of D (Area_D)
# Area_D = Area_adj + 3 * A_opp
# Area_D = 3*s^2 + s^2 * (pi - 9/2 + 3*sqrt(3)/2)
# Area_D = s^2 * (3 + pi - 9/2 + 3*sqrt(3)/2)
# Area_D = s^2 * (pi - 3/2 + 3*sqrt(3)/2)
# To clear the fraction in the bracket, we write:
# Area_D = (s^2 / 2) * (2*pi - 3 + 3*sqrt(3))
print("Step 4: Total Area of Region D")
print("Area(D) = Area_adj + 3 * A_opp")
print("Area(D) = 3*s^2 + s^2 * (pi - 9/2 + 3*sqrt(3)/2) = (s^2 / 2) * (2*pi - 3 + 3*sqrt(3)).")
print("-" * 30)

# Step 5: Final Ratio
# Ratio = Area_D / Area_S
# Ratio = [ (s^2 / 2) * (2*pi - 3 + 3*sqrt(3)) ] / (6 * s^2)
# The s^2 terms cancel out.
# Ratio = (2*pi - 3 + 3*sqrt(3)) / 12
print("Step 5: Final Ratio")
print("The final ratio is Area(D) / Area(S).")
numerator_str = "2*pi - 3 + 3*sqrt(3)"
denominator_str = "12"
final_answer_str = f"({numerator_str}) / {denominator_str}"
print(f"The numerator of the final fraction is: {numerator_str}")
print(f"The denominator of the final fraction is: {denominator_str}")
print(f"\nThe final ratio in its exact form is:\n{final_answer_str}")
