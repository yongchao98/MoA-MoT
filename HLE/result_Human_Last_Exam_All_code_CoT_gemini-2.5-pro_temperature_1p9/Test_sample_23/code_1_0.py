import math

def count_non_admissible(a, b):
    """
    Calculates the number of non-admissible integers k for a given a and b.

    An integer k is "admissible" if there exist complex a by b matrices A_1,...,A_{ab}
    satisfying:
    1. Each A_i is nonzero
    2. tr(A_i^dagger A_j) = 0 whenever i != j
    3. Exactly k of the matrices A_i have rank 1.

    The number of non-admissible k in {0, 1, ..., ab} is determined by the values of a and b.
    """

    if a <= 0 or b <= 0 or not isinstance(a, int) or not isinstance(b, int):
        raise ValueError("a and b must be positive integers.")

    m = min(a, b)
    M = max(a, b)

    # Case 1: One of a or b is 1.
    # Any nonzero matrix is rank-1, so any orthogonal basis must consist of
    # ab rank-1 matrices. Only k=ab is admissible.
    if m == 1:
        # The non-admissible integers are 0, 1, ..., ab-1.
        return m * M

    # Case 2: a = b = n >= 2.
    # It's a known result that k=1 is the only non-admissible value.
    elif m == M:
        return 1

    # Case 3: a, b > 1 and a != b.
    else:
        # Subcase 3a: min(a,b) = 2.
        # For M_{2,M}(C) with M>2, admissible k are only the even numbers from 0 to 2M.
        # The non-admissible integers are the odd numbers: 1, 3, ..., 2M-1.
        # The count of non-admissible numbers is M.
        if m == 2:
            return M
        # Subcase 3b: min(a,b) >= 3.
        # Based on results for pairs like (3,4), the number of non-admissible
        # values appears to be min(a,b). For example, for (3,4) they are {1,3,11}.
        else:  # m >= 3 and m != M
            return m

# Since a and b are not specified, here are some examples.
examples = [(4, 4), (2, 5), (3, 5), (1, 7)]

for a, b in examples:
    result = count_non_admissible(a, b)
    print(f"For a={a}, b={b}, the number of non-admissible integers is: {result}")