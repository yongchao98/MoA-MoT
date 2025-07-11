# The curve is z^2 = 2*x^5 + 2*x^3 + 1.
# We perform a change of variables z = 1 + 2*eta.
# (1 + 2*eta)^2 = 1 + 4*eta + 4*eta^2 = 2*x^5 + 2*x^3 + 1
# This simplifies to 2*eta^2 + 2*eta = x^5 + x^3.
# The special fiber (mod 2) is x^3 * (x+1)^2 = 0.
# In the projective plane, this is X^3 * (X+W)^2 = 0, with components C1: X=0 (mult 3) and C2: X+W=0 (mult 2).
# They intersect at P=(0:1:0).
# The local equation of the model at P is F(X,W) = X^5 + X^3*W^2 - 2*W^3 - 2*W^4 = 0.
# The term of lowest total degree in F(X,W) is -2*W^3.
# The degree of this term is 3.
# The 2-adic valuation of its coefficient is v(-2) = 1.
# A formula for the thickness of the resulting singularity in the stable reduction is given by the degree of this term minus 1.
degree_of_lowest_term = 3
valuation_of_coefficient = 1

# The thickness is the degree of the lowest-degree part of the equation minus one.
thickness = degree_of_lowest_term - valuation_of_coefficient

# Output the equation
print(f"The degree of the lowest order term is {degree_of_lowest_term}.")
print(f"The valuation of the coefficient of this term is {valuation_of_coefficient}.")
print(f"The thickness is calculated as: {degree_of_lowest_term} - {valuation_of_coefficient} = {thickness}")
