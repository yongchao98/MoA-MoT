from fractions import Fraction

def solve():
    """
    This function calculates the final sum based on the identified set S.
    The set S consists of three disjoint families of pairs (i, j):
    1. (k, k) for k >= 1
    2. (k, 2k) for k >= 1
    3. (2k, k) for k >= 1
    The total sum is computed by summing the corresponding geometric series for each family.
    """

    # Helper function to calculate the sum of an infinite geometric series
    # Sum = a / (1 - r), where a is the first term and r is the common ratio.
    def geometric_series_sum(a, r):
        return a / (1 - r)

    # Sum for the family S1 = {(k, k) | k >= 1}
    # This corresponds to the series Sum_{k=1 to inf} 1/2^(k+k) = Sum (1/4)^k
    # Here, a = 1/4 and r = 1/4.
    sum1_a = Fraction(1, 4)
    sum1_r = Fraction(1, 4)
    sum1 = geometric_series_sum(sum1_a, sum1_r)

    # Sum for the family S2 = {(k, 2k) | k >= 1}
    # This corresponds to the series Sum_{k=1 to inf} 1/2^(k+2k) = Sum (1/8)^k
    # Here, a = 1/8 and r = 1/8.
    sum2_a = Fraction(1, 8)
    sum2_r = Fraction(1, 8)
    sum2 = geometric_series_sum(sum2_a, sum2_r)

    # Sum for the family S3 = {(2k, k) | k >= 1}
    # This corresponds to the series Sum_{k=1 to inf} 1/2^(2k+k) = Sum (1/8)^k
    # This is identical to the sum for S2.
    sum3 = sum2

    # The total sum is the sum of the three disjoint series.
    total_sum = sum1 + sum2 + sum3

    # Output the details of the calculation.
    print("The set S of pairs (i, j) is the union of three disjoint sets:")
    print("S1 = {(k, k) | k>=1}, S2 = {(k, 2k) | k>=1}, S3 = {(2k, k) | k>=1}")
    print("\nThe total sum is Sum(S1) + Sum(S2) + Sum(S3).")
    print(f"\nSum over S1: Sum_{{k=1 to inf}} 1/2^(k+k) = {sum1}")
    print(f"Sum over S2: Sum_{{k=1 to inf}} 1/2^(k+2k) = {sum2}")
    print(f"Sum over S3: Sum_{{k=1 to inf}} 1/2^(2k+k) = {sum3}")
    
    print("\nFinal sum equation:")
    print(f"{sum1} + {sum2} + {sum3} = {total_sum}")

solve()