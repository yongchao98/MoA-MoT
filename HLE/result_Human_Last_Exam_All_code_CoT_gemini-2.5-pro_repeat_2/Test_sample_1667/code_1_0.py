import fractions

# This script calculates the probability based on the analytical solution derived above.

# The final probability 'p' is given by the formula:
# p = (2 - P_tri) / 4
# where P_tri is the probability that the convex hull of 4 random points in a unit square is a triangle.

# P_tri is calculated as 4 * E_area, where E_area is the expected area of a
# random triangle in a unit square. This is a known result from geometric probability.
E_area_num = 11
E_area_den = 144
E_area = fractions.Fraction(E_area_num, E_area_den)

# The equation for P_tri uses the number 4 and the expected area.
P_tri_factor = 4
P_tri = P_tri_factor * E_area

# The final equation for the probability 'p' uses the numbers 2 and 4, and the value of P_tri.
p_numerator_term = 2
p_denominator_term = 4
p = (p_numerator_term - P_tri) / p_denominator_term

# Output the steps of the calculation as a clear narrative.
print(f"The formula for the probability is p = (2 - P_tri) / 4")
print(f"P_tri is the probability that 4 random points form a triangle.")
print(f"P_tri is calculated as {P_tri_factor} * (Expected Area of a random triangle).")
print(f"The expected area of a random triangle in a unit square is a known value: {E_area_num}/{E_area_den}.")
print(f"So, P_tri = {P_tri_factor} * {E_area_num}/{E_area_den} = {P_tri.numerator}/{P_tri.denominator}.")
print(f"Substituting this into the formula for p:")
print(f"p = ({p_numerator_term} - {P_tri.numerator}/{P_tri.denominator}) / {p_denominator_term}")
print(f"p = (({p_numerator_term*P_tri.denominator - P_tri.numerator})/{P_tri.denominator}) / {p_denominator_term}")
print(f"p = {p.numerator}/{p.denominator}")
print(f"\nThe final probability is {p.numerator}/{p.denominator}, which is approximately {float(p):.4f}.")
