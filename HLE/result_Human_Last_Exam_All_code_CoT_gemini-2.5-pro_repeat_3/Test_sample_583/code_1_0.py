import math

# The problem is to evaluate the limit L = lim_{n->inf} f(n) / (n * log2(n)).
#
# From our analysis, we determined that the function f(n) has an upper bound
# that grows linearly with n, i.e., f(n) = O(n).
#
# The limit can be written as:
# L = lim_{n->inf} O(n) / (n * log2(n))
#
# This simplifies to:
# L = lim_{n->inf} C / log2(n) for some constant C.
#
# As n approaches infinity, log2(n) also approaches infinity.
# Therefore, the limit is 0.

limit_value = 0
print(limit_value)