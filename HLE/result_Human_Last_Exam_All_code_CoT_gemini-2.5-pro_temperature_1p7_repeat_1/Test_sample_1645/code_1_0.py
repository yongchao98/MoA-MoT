# The problem asks for the smallest non-negative integer n such that the property (Rn)
# is not preserved by completion of a noetherian local ring.

# (Rn) means "regular in codimension <= n". A ring A has (Rn) if for every
# prime ideal p of A with height(p) <= n, the localization A_p is a regular local ring.

# For n = 0: The property (R0) is always preserved under m-adic completion for a
# noetherian local ring. (R0) means that the localization at every minimal prime is a field.
# Thus, n cannot be 0.

# For n = 1: The property (R1) is not always preserved. There exists a famous
# counterexample, constructed by Nagata, of a 2-dimensional regular local ring A
# whose completion Â is not normal.
# - Since A is a regular local ring, it satisfies (Rn) for all n, in particular for n=1.
# - A ring is normal if and only if it satisfies properties (R1) and (S2) (Serre's Criterion).
# - Since A is a regular local ring, it is Cohen-Macaulay. The completion Â is also
#   Cohen-Macaulay, which implies Â satisfies (S2).
# - Since Â is not normal but satisfies (S2), it must fail to satisfy (R1).
# This means there is a prime ideal q in Â of height at most 1 such that (Â)_q
# is not regular.

# So, we have an example where A satisfies (R1), but Â does not.

# Therefore, the smallest non-negative integer n for which (Rn) is not
# preserved by completion is 1.

n = 1

# There is no equation, so we will just print the final answer for n.
print(n)