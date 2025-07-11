import math

def solve_pm(m):
    """
    Calculates the probability P_m for a given positive integer m.

    m is a positive integer, and the sequence a_1, ..., a_{4m+2} is an arithmetic sequence.
    P_m is the probability that the sequence is (i,j)-divisible.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: m must be a positive integer.")
        return

    # The number of pairs (i,j) that make the sequence (i,j)-divisible
    # is 3 if m is odd, and 5 if m is even.
    # This can be written as 4 + (-1)^m.
    numerator = 4 + (-1)**m

    # The total number of ways to choose i and j from {1, ..., 4m+2} with i < j
    # is C(4m+2, 2).
    denominator = (2 * m + 1) * (4 * m + 1)

    # Use math.gcd to simplify the fraction
    common_divisor = math.gcd(numerator, denominator)
    simplified_num = numerator // common_divisor
    simplified_den = denominator // common_divisor

    print(f"For m = {m}:")
    print(f"The number of valid (i,j) pairs is 4 + (-1)^{m} = {numerator}")
    print(f"The total number of (i,j) pairs is (2m+1)(4m+1) = (2*{m}+1)(4*{m}+1) = {denominator}")
    print(f"The probability P_m = {numerator} / {denominator}")
    print(f"Simplified P_m = {simplified_num} / {simplified_den}")

# Example usage of the function. You can change the value of m.
# Let's calculate for m=1 (an odd number)
solve_pm(1)
print("-" * 20)
# Let's calculate for m=2 (an even number)
solve_pm(2)
print("-" * 20)
# Let's calculate for m=3 (an odd number)
solve_pm(3)
