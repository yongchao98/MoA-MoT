from fractions import Fraction

def calculate_sum():
    """
    This function calculates the required sum based on the derivation.
    The set S consists of pairs (i, j) where j=i, j=2i, or i=2j.
    The total sum is computed by summing over these three disjoint sets of pairs.
    """

    # The sum of an infinite geometric series sum_{k=1 to inf} r^k is r / (1-r).
    def geometric_series_sum(r_num, r_den):
        r = Fraction(r_num, r_den)
        return r / (1 - r)

    # Case 1: Pairs (k, k) for k >= 1. The sum is Sum_{k=1 to inf} (1/2^(2k)) = Sum (1/4)^k.
    # The ratio is r = 1/4.
    sum_case1 = geometric_series_sum(1, 4)

    # Case 2: Pairs (k, 2k) for k >= 1. The sum is Sum_{k=1 to inf} (1/2^(3k)) = Sum (1/8)^k.
    # The ratio is r = 1/8.
    sum_case2 = geometric_series_sum(1, 8)
    
    # Case 3: Pairs (2k, k) for k >= 1. The sum is Sum_{k=1 to inf} (1/2^(3k)) = Sum (1/8)^k.
    # This is identical to case 2.
    sum_case3 = sum_case2

    # The total sum is the sum of the three cases.
    total_sum = sum_case1 + sum_case2 + sum_case3

    # Print the result in the required format.
    print("The sum is composed of three parts based on the pairs (i, j) in S:")
    print(f"1. For pairs where j=i: Sum = {sum_case1}")
    print(f"2. For pairs where j=2i: Sum = {sum_case2}")
    print(f"3. For pairs where i=2j: Sum = {sum_case3}")
    print("\nThe final equation is:")
    print(f"{sum_case1} + {sum_case2} + {sum_case3} = {total_sum}")

calculate_sum()