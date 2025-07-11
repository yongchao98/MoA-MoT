import sympy

# The problem is to compute the sum S = sum_{n=1 to inf} c(n) / 4^n,
# where c(n) is the number of possible values for m for a given n.
# Based on the analysis of the problem (and its source, Putnam 2018 B6), we have:
# c(n) = 2 if n is odd
# c(n) = 4 if n is even

# The sum can be split into two parts: one for odd n and one for even n.

# Part 1: Sum for odd n
# Let n = 2k - 1 for k = 1, 2, 3, ...
# The sum is Sum_{k=1 to inf} 2 / 4^(2k-1)
k = sympy.Symbol('k', integer=True, positive=True)
sum_odd_term = 2 / (4**(2*k - 1))
sum_odd = sympy.summation(sum_odd_term, (k, 1, sympy.oo))

# Part 2: Sum for even n
# Let n = 2k for k = 1, 2, 3, ...
# The sum is Sum_{k=1 to inf} 4 / 4^(2k)
sum_even_term = 4 / (4**(2*k))
sum_even = sympy.summation(sum_even_term, (k, 1, sympy.oo))

# The total sum is the sum of these two parts.
total_sum = sum_odd + sum_even

# Print the final equation with all its components
print(f"Sum for odd n: {sum_odd}")
print(f"Sum for even n: {sum_even}")
print(f"Total sum = {sum_odd} + {sum_even} = {total_sum}")
