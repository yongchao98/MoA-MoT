import sys

def solve():
    """
    Calculates the number of non-admissible integers k for a x b matrices.

    An integer k is "admissible" if there exist ab non-zero complex a x b
    matrices A_1, ..., A_{ab} satisfying the following conditions:
    1. Each A_i is nonzero
    2. tr(A_i^dagger A_j) = 0 whenever i != j
    3. Exactly k of the matrices A_i have rank 1.

    The number of integers in the range 0, 1, ..., ab that are not
    admissible is calculated.
    """
    # You can change the positive integer values for a and b here.
    a = 4
    b = 6

    if not isinstance(a, int) or not isinstance(b, int) or a <= 0 or b <= 0:
        print("Error: Please provide positive integers for a and b.")
        return

    # Let m = min(a,b) and M = max(a,b).
    m = min(a, b)
    M = max(a, b)

    num_non_admissible = 0

    if m == 1:
        # Case 1: If min(a,b) = 1, any non-zero matrix has rank 1.
        # An orthogonal basis must consist of non-zero matrices, so all ab
        # matrices have rank 1. k = ab is the only admissible value.
        # The non-admissible values are {0, 1, ..., ab-1}.
        num_non_admissible = m * M
        print_explanation = f"min(a,b) is 1. The number of non-admissible integers is ab = {m}*{M}."
    elif m == M:
        # Case 2: If a = b >= 2, all values of k are admissible.
        num_non_admissible = 0
        print_explanation = f"a=b>=2. The number of non-admissible integers is 0."
    else:  # This covers the case 2 <= m < M
        # Case 3: If 2 <= min(a,b) < max(a,b), the number of non-admissible
        # integers is min(a,b) - 1.
        num_non_admissible = m - 1
        print_explanation = f"2<=min(a,b)<max(a,b). The number is min(a,b) - 1 = {m}-1."

    print(f"Given a = {a} and b = {b}:")
    print(print_explanation)
    print(f"Result: {num_non_admissible}")


solve()