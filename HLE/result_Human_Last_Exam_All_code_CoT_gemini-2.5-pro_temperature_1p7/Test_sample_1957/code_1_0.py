# Plan:
# The problem asks for min(1000*m + n).
# My analysis shows that the set S is the entire space of 5-point tuples.
# This means det(F) must be the zero polynomial.
# We need to find the smallest positive integers m and n for which
# such an F can be constructed.

# The smallest possible dimension n for a square matrix F to have a
# determinant that is identically zero (while its entries are non-zero polynomials)
# is n=2. For example, a matrix with two identical rows:
# F = [[p, q], [p, q]], where p and q are polynomials.
# det(F) = p*q - q*p = 0.

# The degree 'm' of the polynomial map F is the maximum degree among its
# polynomial entries. The problem states m must be a positive integer, so the
# minimum value for m is 1.

# We can choose p and q to be simple polynomials of degree 1, for example,
# the coordinates of one of the points.

m = 1  # Smallest positive integer for the degree.
n = 2  # Smallest dimension for det(F) to be identically zero.

# Calculate the result.
result = 1000 * m + n

# The problem asks to output the numbers in the final equation.
print(f"The minimum value is obtained with m = {m} and n = {n}.")
print(f"The calculation is: 1000 * {m} + {n} = {result}")
