import math

def solve():
    """
    Calculates the number of non-admissible integers for given positive integers a and b.
    """
    # Let's define the positive integers a and b.
    # You can change these values to test other cases.
    a = 4
    b = 5

    print(f"For a = {a} and b = {b}:")

    # The total number of matrices is ab.
    prod = a * b

    # Case 1: If a or b is 1, the matrices are vectors.
    # Any non-zero vector has rank 1. Since all A_i must be non-zero,
    # all ab matrices in the basis must have rank 1.
    # So, k must be ab.
    if min(a, b) == 1:
        result = prod
        # The non-admissible k are 0, 1, ..., ab-1. There are ab of them.
        print(f"The number of non-admissible integers is {a} * {b} = {result}")
    else:
        # Case 2: a,b >= 2.
        # An integer k is admissible if and only if it has the same parity as ab.
        # The non-admissible integers are those with the opposite parity.
        if prod % 2 == 0:
            # If ab is even, non-admissible integers are the odd ones.
            # The number of odd integers in [0, 1, ..., ab] is ab/2.
            result = prod // 2
            print(f"ab is even. The number of non-admissible integers is {a} * {b} / 2 = {result}")
        else:
            # If ab is odd, non-admissible integers are the even ones.
            # The number of even integers in [0, 1, ..., ab] is (ab+1)/2.
            result = (prod + 1) // 2
            print(f"ab is odd. The number of non-admissible integers is ({a} * {b} + 1) / 2 = {result}")

solve()
<<<10>>>