# The user is asking for the smallest non-negative integer n such that the property (Rn)
# is not necessarily preserved when taking the completion of a Noetherian local ring.
#
# (Rn) is a property of a Noetherian ring A, which states that for any prime ideal p
# in A with height(p) <= n, the localization A_p is a regular local ring.
#
# The question is: What is the minimum n >= 0 for which there exists a ring A
# satisfying (Rn), but its completion Â does not?

# Let's analyze the problem for small values of n.

# Case n = 0:
# We check if (R0) is preserved under completion.
# Assume A has (R0). Let P be a prime ideal in Â with height(P) = 0.
# Let p = P intersect A. Since the completion map A -> Â is flat, we have height(p) <= height(P), so height(p) = 0.
# Since A has (R0), A_p is a regular local ring. A regular local ring of dimension 0 is a field.
# The formal fibers of a field are regular.
# The regularity of Â_P is equivalent to the regularity of A_p and the regularity of the formal fiber.
# Since both are regular, Â_P is regular.
# This holds for all primes P of height 0, so Â has (R0).
# Thus, (R0) is always preserved. So n is not 0.

# Case n = 1:
# We check if (R1) is preserved under completion.
# Assume A has (R1). Let P be a prime ideal in Â with height(P) <= 1.
# Let p = P intersect A. Due to flatness, height(p) <= height(P) <= 1.
# Since A has (R1), A_p is a regular local ring of dimension at most 1.
# A regular local ring of dimension 0 (a field) or 1 (a DVR) is always a Nagata ring,
# and therefore a G-ring. This means its formal fibers are geometrically regular.
# Since A_p is regular and its formal fibers are regular, Â_P must be regular.
# This holds for all primes P of height at most 1, so Â has (R1).
# Thus, (R1) is always preserved. So n is not 1.

# Case n = 2:
# For n=2, the situation changes. A regular local ring of dimension 2 is not necessarily
# a Nagata ring (and thus not a G-ring). This means its formal fibers are not guaranteed
# to be regular.
# Indeed, counterexamples for n=2 have been constructed by prominent algebraists.
# Nagata, in the 1950s, and later Ogoma (1980), constructed examples of Noetherian local rings A
# which satisfy (R2), but whose completions Â do not satisfy (R2).
# These constructions are highly non-trivial.
# For example, one can construct a 3-dimensional normal (hence R1) Noetherian local domain A
# which is also (R2), but its completion Â has a non-regular prime of height 2.
# This shows that (R2) is not preserved.

# Conclusion:
# Since (R0) and (R1) are always preserved under completion, and a counterexample exists for (R2),
# the smallest such non-negative integer n is 2.

smallest_n = 2
print(smallest_n)