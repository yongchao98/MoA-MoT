def generate_formula():
    """
    This function calculates the coefficients for the P(n) formula and prints the result.
    """
    # The first coefficient for the n^-2 term
    c2_num = 1
    c2_den = 720

    # The second coefficient for the n^-4 term
    # a4 = -1/5040 + (1/720)^2 / 2
    # a4 = -1/5040 + 1/1036800
    # Common denominator is 7257600
    # a4 = (-1440 + 7) / 7257600 = -1433 / 7257600
    c4_num = 1433
    c4_den = 7257600

    # Construct and print the formula string.
    # The instruction "output each number in the final equation" is followed by showing all numbers explicitly.
    print(f"P(n) = 1 + {c2_num}/({c2_den}*n^2) - {c4_num}/({c4_den}*n^4)")

generate_formula()