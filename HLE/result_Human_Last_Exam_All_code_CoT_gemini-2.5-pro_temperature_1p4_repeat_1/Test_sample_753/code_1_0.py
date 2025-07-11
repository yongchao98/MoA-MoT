import math
from decimal import Decimal, getcontext

def calculate_cardinality_digits():
    """
    This function calculates the first 40 digits of the cardinality of Theta^{-1}(lambda) for m=3.
    """
    # Set precision for decimal calculations. A high precision is needed for the logarithm of a large number.
    getcontext().prec = 100

    # Step 1: Determine the value of n for m=3.
    m = 3
    n = sum(k * (m + 1 - k) for k in range(1, m + 1))
    
    # Step 2: Determine the partition lambda and calculate the size of its conjugacy class |C_lambda|.
    # For m=3, the partition lambda is (3^1, 2^2, 1^3).
    # This corresponds to a permutation in Sigma_10 with cycle type (1,1,1,2,2,3).
    # The size of the conjugacy class |C_lambda| in Sigma_n is n! / (product of k^(n_k) * n_k!).
    
    # Calculate n!
    n_factorial = math.factorial(n)
    
    # The denominator for |C_lambda| calculation for lambda = (3^1, 2^2, 1^3) is:
    # (3^1 * 1!) * (2^2 * 2!) * (1^3 * 3!) = 3 * 8 * 6 = 144
    denominator = (3**1 * math.factorial(1)) * (2**2 * math.factorial(2)) * (1**3 * math.factorial(3))
    
    # Calculate the size of the conjugacy class.
    C_lambda_size = n_factorial // denominator
    
    # The final formula for the cardinality is |C_lambda| * (n!)^(n!).
    # We print the components of this formula as requested.
    print("The calculation is based on the formula: |C_lambda| * (n!)^(n!)")
    print(f"For m = 3, we have n = {n}.")
    print(f"n! = {n_factorial}")
    print(f"|C_lambda| = {C_lambda_size}")
    print(f"The final number is given by the expression: {C_lambda_size} * ({n_factorial})^({n_factorial})")
    
    # Step 3: Calculate the first 40 digits of the result using logarithms.
    # Let N be the result. We calculate log10(N) = log10(|C_lambda|) + n! * log10(n!).
    
    val_C = Decimal(C_lambda_size)
    val_n_fact = Decimal(n_factorial)
    
    log10_N = val_C.log10() + val_n_fact * val_n_fact.log10()
    
    # The first digits of N are determined by 10^f, where f is the fractional part of log10(N).
    # For a number N = 10^(M+f) = (10^f) * 10^M, 10^f gives the significant digits.
    f = log10_N - log10_N.to_integral_value(rounding=math.floor)
    
    first_digits_val = Decimal(10)**f
    
    # Format the resulting Decimal object into a string of the first 40 digits.
    s = str(first_digits_val)
    first_40_digits = (s[0] + s[2:])[:40]
    
    print("\nThe first 40 digits of the cardinality are:")
    print(first_40_digits)

calculate_cardinality_digits()