# This problem is a theoretical question from the field of topology.
# The solution relies on a known mathematical theorem regarding the Stone-Cech compactification.
#
# Let f: R -> R be a continuous function and F be its extension to the Stone-Cech
# compactification of R, denoted beta(R).
# The Stone-Cech remainder is R* = beta(R) \ R.
# We are looking for the smallest possible nonzero number of fixed points of F in R*.
#
# A theorem by A. Szymański (1994) states that the number of fixed points of F
# in R* can only be 0, 1, 2, or 2^c (where c is the cardinality of the continuum).
#
# The question asks for the smallest *nonzero* number among these possibilities.
# The set of nonzero possibilities is {1, 2, 2^c}.
#
# The smallest value in this set is 1.
#
# It has been shown that a function f can be constructed such that its extension F
# has exactly one fixed point in R*. Therefore, the number 1 is an achievable value.
#
# Thus, the smallest possible nonzero number of fixed points is 1.

smallest_nonzero_fixed_points = 1

# The final equation is simply the value of the answer.
# Outputting the number in the final equation:
print(smallest_nonzero_fixed_points)