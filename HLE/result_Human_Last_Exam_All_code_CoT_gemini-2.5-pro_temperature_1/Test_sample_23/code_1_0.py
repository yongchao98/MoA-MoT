import math

def solve():
    """
    Calculates the number of non-admissible integers k for given a and b.

    An integer k is "admissible" if there exist complex a by b matrices A_1,...,A_{ab}
    satisfying the following conditions:
    1. Each A_i is nonzero
    2. tr(A_i^dagger A_j) = 0 whenever i != j
    3. exactly k of the matrices A_i have rank 1.

    The number of non-admissible integers in the range [0, 1, ..., ab] depends on a and b.
    """

    # Let a,b be positive integers.
    # You can change these values to test with other numbers.
    a = 3
    b = 5

    product = a * b

    # Case 1: If either dimension is 1, any non-zero matrix has rank 1.
    # An orthogonal basis must consist of ab non-zero matrices, all of which will have rank 1.
    # So, k must be ab. All other values are not admissible.
    if min(a, b) == 1:
        num_not_admissible = product
        not_admissible_list = list(range(product))
    # Case 2: Both dimensions are 2 or more.
    # It can be shown that the set of admissible k are all integers in [0, ab]
    # that have the same parity as ab.
    else:
        # If ab is even, admissible k are {0, 2, ..., ab}.
        # Not admissible k are {1, 3, ..., ab-1}.
        if product % 2 == 0:
            num_not_admissible = product // 2
            not_admissible_list = [i for i in range(1, product, 2)]
        # If ab is odd, admissible k are {1, 3, ..., ab}.
        # Not admissible k are {0, 2, ..., ab-1}.
        else:
            num_not_admissible = (product + 1) // 2
            not_admissible_list = [i for i in range(0, product, 2)]

    print(f"For a={a}, b={b}:")
    print(f"The non-admissible integers are: {not_admissible_list}")
    print(f"The number of non-admissible integers is: {num_not_admissible}")

solve()