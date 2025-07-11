# Plan:
# The geometric condition described in the problem is satisfied by any set of 5 points (A, B, C, D, X).
# This means the set S is the entire space of possible inputs.
# The condition `det(F) = 0` must hold for all inputs, which implies det(F) must be the zero polynomial.
# We need to find the smallest positive integers n and m for which we can construct
# an n x n matrix F with polynomial entries of degree at most m such that det(F) is identically zero.

# n must be at least 2 for the determinant of a non-zero matrix of polynomials to be the zero polynomial.
# For n=2, we can have a matrix like F = [[P, Q], [P, Q]], where P, Q are polynomials.
# The determinant is then P*Q - Q*P = 0.
# The minimum positive integer value for m (the degree of the polynomial entries) is 1.
# We can choose P and Q to be polynomials of degree 1 (e.g., coordinates of the points).

# So, the minimal positive integers are m=1 and n=2.

m = 1
n = 2

# The expression to minimize is 1000m + n.
result = 1000 * m + n

# The final equation as requested by the prompt format
print(f"m = {m}")
print(f"n = {n}")
print(f"1000 * {m} + {n} = {result}")
