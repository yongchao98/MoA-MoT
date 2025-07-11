import math

def get_smallest_prime_factor(n):
    """
    Returns the smallest prime factor of a given integer n.
    """
    if n % 2 == 0:
        return 2
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        if n % i == 0:
            return i
    return n

def solve():
    """
    Calculates the number of positive integers n <= lcm(1, ..., 100)
    that have distinct remainders when divided by k = 2, ..., 100.

    It is a known result that such integers must be of the form
    n = L - c, where L = lcm(1, ..., 100).
    This function counts the number of valid values for c.
    A value c is valid if for all composite k from 4 to 100,
    the remainder (-c) mod k is greater than or equal to the largest
    proper divisor of k.
    """
    solution_count = 0
    # We only need to check for small values of c. The conditions for c
    # are periodic, and a failure is likely to occur for a small k.
    # Checking up to 100 is sufficient to establish a pattern.
    for c in range(1, 101):
        is_valid_c = True
        for k in range(4, 101):
            # Check if k is composite
            p_min = get_smallest_prime_factor(k)
            if p_min < k:  # k is composite
                largest_proper_divisor = k // p_min
                remainder = (-c) % k
                if remainder < largest_proper_divisor:
                    is_valid_c = False
                    break
        
        if is_valid_c:
            solution_count += 1
            # print(f"Found a solution for c = {c}")

    print(f"The number of such positive integers is {solution_count}.")

solve()
