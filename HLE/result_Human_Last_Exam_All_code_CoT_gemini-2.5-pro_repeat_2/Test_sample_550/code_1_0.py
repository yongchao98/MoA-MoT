import math

# This script calculates the dimension of the ninth cohomology group H^9(M, Q).
# The space M is the complement of a specific quaternionic hyperplane arrangement in H^4.
# The calculation relies on identifying the arrangement and using known results from the theory of free arrangements.

# The key insight is that the given arrangement of 36 hyperplanes is a "free arrangement"
# associated with the complex reflection group G_31. The properties of this arrangement
# are encoded in its "exponents".
exponents = [3, 7, 11, 15]
d1, d2, d3, d4 = exponents

# For a free quaternionic arrangement, the Betti numbers b_k = dim H^k(M, Q) are the
# coefficients of the Poincaré polynomial P(t) = (1 + d1*t^3)(1 + d2*t^3)(1 + d3*t^3)(1 + d4*t^3).
# We are looking for b_9, the coefficient of t^9.
# This coefficient is the third elementary symmetric polynomial of the exponents.
# b_9 = e_3(d1, d2, d3, d4) = d1*d2*d3 + d1*d2*d4 + d1*d3*d4 + d2*d3*d4

# Calculate each term of the sum
term1 = d1 * d2 * d3
term2 = d1 * d2 * d4
term3 = d1 * d3 * d4
term4 = d2 * d3 * d4

# Calculate the final result
result = term1 + term2 + term3 + term4

print(f"The dimension of H^9(M, Q) is given by the third elementary symmetric polynomial of the exponents {exponents}.")
print(f"The calculation is: {d1}*{d2}*{d3} + {d1}*{d2}*{d4} + {d1}*{d3}*{d4} + {d2}*{d3}*{d4}")
print(f"The values of the terms are: {term1} + {term2} + {term3} + {term4}")
print(f"The final dimension is: {result}")