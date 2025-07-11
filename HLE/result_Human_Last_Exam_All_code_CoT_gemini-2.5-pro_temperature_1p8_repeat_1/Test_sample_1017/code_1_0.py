from fractions import Fraction

# Define the translation numbers for g and h
tau_g = Fraction(2, 27)
tau_h = Fraction(16, 27)

# The stable commutator length (scl) of an element f in a group like Homeo+(S^1)
# is related to the norm of its class in the abelianization H_1(G),
# which is given by the absolute value of its translation number tau(f).
# scl(f) = 1/2 * ||[f]|| = 1/2 * |tau(f)|

# For a free product G1 * G2, the scl norm on the abelianization H_1(G1) x H_1(G2)
# is the L^1 norm. For an element g1*h2, the calculation is:
# scl(g1*h2) = 1/2 * ||([g1], [h2])|| = 1/2 * (||[g1]|| + ||[h2]||)
#            = 1/2 * (|tau(g1)| + |tau(h2)|)

# Perform the calculation
scl_result = (abs(tau_g) + abs(tau_h)) / 2

# Print the final equation with all numbers
print(f"scl(g1*h2) = (1/2) * (|{tau_g.numerator}/{tau_g.denominator}| + |{tau_h.numerator}/{tau_h.denominator}|) = {scl_result.numerator}/{scl_result.denominator}")