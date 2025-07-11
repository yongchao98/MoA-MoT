# The problem asks for the smallest non-negative integer n such that the property (Rn)
# is not preserved by completion of a Noetherian local ring.
#
# (Rn) means "regular in codimension at most n".
#
# Based on the theory of commutative algebra:
# - For n = 0, the property (R0) is always preserved by completion. This can be proven using
#   the flatness of the completion map.
# - For n = 1, the property (R1) is NOT always preserved. Famous counterexamples
#   were constructed by Nagata, who found a normal (and thus R1) local domain A
#   whose completion Â is not (R1).
#
# Therefore, the smallest such non-negative integer is 1.

smallest_n = 1

# There is no equation, so we will just print the final answer.
print(smallest_n)
