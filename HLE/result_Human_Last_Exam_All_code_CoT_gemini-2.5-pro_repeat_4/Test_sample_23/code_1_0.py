import math

def count_non_admissible_integers():
    """
    Calculates the number of non-admissible integers k for given matrix dimensions a and b.

    An integer k is "admissible" if there exists an orthogonal basis for the space
    of a x b complex matrices where exactly k of the basis matrices have rank 1.
    The function determines the number of integers in the range {0, 1, ..., ab}
    that are NOT admissible.
    """
    # Let a, b be positive integers.
    # You can change these values to the desired dimensions.
    a = 3
    b = 5

    print(f"For matrices of size a={a}, b={b}:")

    # Case 1: At least one dimension is 1 (matrices are essentially vectors).
    if a == 1 or b == 1:
        # In this case, any non-zero matrix has rank 1. An orthogonal basis must consist
        # of non-zero matrices, so all 'ab' basis matrices must have rank 1.
        # The only admissible value for k is 'ab'.
        # The number of non-admissible integers is (ab + 1) - 1 = ab = max(a, b).
        result = max(a, b)

        print("This falls into Case 1 (min(a, b) = 1).")
        print("The number of non-admissible integers is max(a, b).")
        # Output the equation with the specific numbers
        print(f"max({a}, {b}) = {result}")

    # Case 2: Both dimensions are 2 or greater.
    else:
        # In this case, the set of admissible k's are all integers from a minimum
        # value up to 'ab' that have the same parity as 'ab'. The minimum is 'ab mod 2'.
        # The number of non-admissible integers simplifies to ceil(ab / 2).
        ab = a * b
        # In integer arithmetic, ceil(n / 2) is calculated as (n + 1) // 2.
        result = (ab + 1) // 2

        print("This falls into Case 2 (min(a, b) >= 2).")
        print("The number of non-admissible integers is ceil(a*b / 2).")
        # Output the equation with the specific numbers
        print(f"({a} * {b} + 1) // 2 = ({ab} + 1) // 2 = {result}")

# Run the function to get the answer.
count_non_admissible_integers()