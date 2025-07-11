from fractions import Fraction

def solve_sum():
    """
    This function calculates the specified sum based on the derived conditions for the set S.

    The set S consists of pairs (i, j) of positive integers that fall into three disjoint categories:
    1. (k, k) for any positive integer k.
    2. (2k, k) for any positive integer k.
    3. (k, 2k) for any positive integer k.

    The total sum is the sum of three geometric series corresponding to these categories.
    """

    # Sum for pairs (k, k): sum_{k=1 to inf} 1/2^(2k) = sum (1/4)^k
    # Geometric series with a=1/4, r=1/4. Sum = a / (1-r)
    sum1_val = Fraction(1, 4) / (1 - Fraction(1, 4))

    # Sum for pairs (2k, k): sum_{k=1 to inf} 1/2^(3k) = sum (1/8)^k
    # Geometric series with a=1/8, r=1/8. Sum = a / (1-r)
    sum2_val = Fraction(1, 8) / (1 - Fraction(1, 8))

    # Sum for pairs (k, 2k) is the same as for (2k, k)
    sum3_val = sum2_val

    # The total sum is the sum of the three parts
    total_sum = sum1_val + sum2_val + sum3_val

    # Print the final equation with all its components as requested
    print(f"Sum over pairs (k,k) for k>=1: {sum1_val}")
    print(f"Sum over pairs (2k,k) for k>=1: {sum2_val}")
    print(f"Sum over pairs (k,2k) for k>=1: {sum3_val}")
    print(f"Final Equation: {sum1_val} + {sum2_val} + {sum3_val} = {total_sum}")

if __name__ == '__main__':
    solve_sum()