from fractions import Fraction

def solve():
    """
    Calculates the sum based on the three disjoint sets of pairs (i, j).
    """

    # Sum for pairs of the form (k, k)
    # This is the sum of the geometric series Sum_{k=1 to inf} (1/4)^k
    # First term a = 1/4, common ratio r = 1/4
    # Sum = a / (1 - r)
    sum_set1 = Fraction(1, 4) / (1 - Fraction(1, 4))

    # Sum for pairs of the form (k, 2k) and (2k, k)
    # This is the sum of the geometric series Sum_{k=1 to inf} (1/8)^k
    # First term a = 1/8, common ratio r = 1/8
    # Sum = a / (1 - r)
    sum_set2 = Fraction(1, 8) / (1 - Fraction(1, 8))
    sum_set3 = sum_set2

    # The total sum is the sum of the sums from the three disjoint sets.
    total_sum = sum_set1 + sum_set2 + sum_set3
    
    # Print the final equation with each number.
    print(f"{sum_set1} + {sum_set2} + {sum_set3} = {total_sum}")

solve()